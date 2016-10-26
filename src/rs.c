#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include "fec.h"

#include "fec.c"

#define DATA_SHARDS_MAX (255)
#define RS_MALLOC(x)    malloc(x)
#define RS_FREE(x)      free(x)
#define RS_CALLOC(n, x) calloc(n, x)

typedef struct _reed_solomon {
    int data_shards;
    int parity_shards;
    int shards;
    gf* m;
    gf* parity;
} reed_solomon;

/* only for test */
static void print_buf(gf* buf, char *fmt, size_t len) {
    size_t i = 0;
    while(i < len) {
        printf(fmt, buf[i]);
        i++;
        if((i % 16) == 0) {
            printf("\n");
        }
    }
    printf("\n");
}

static void print_int(int* buf, char *fmt, size_t len) {
    size_t i = 0;
    while(i < len) {
        printf(fmt, buf[i]);
        i++;
        if((i % 16) == 0) {
            printf("\n");
        }
    }
    printf("\n");
}

static void print_matrix1(gf* matrix, int nrows, int ncols) {
    int i, j;
    printf("matrix (%d,%d):\n", nrows, ncols);
    for(i = 0; i < nrows; i++) {
        for(j = 0; j < ncols; j++) {
            printf("%6d ", matrix[i*ncols + j]);
        }
        printf("\n");
    }
}

static void print_matrix2(gf** matrix, int nrows, int ncols) {
    int i, j;
    printf("matrix (%d,%d):\n", nrows, ncols);
    for(i = 0; i < nrows; i++) {
        for(j = 0; j < ncols; j++) {
            printf("%6d ", matrix[i][j]);
        }
        printf("\n");
    }
}

/* y = a**n */
static gf galExp(gf a, gf n) {
    int logA;
    int logResult;
    if(0 == n) {
        return 1;
    }
    if(0 == a) {
        return 0;
    }
    logA = gf_log[a];
    logResult = logA * n;
    while(logResult >= 255) {
        logResult -= 255;
    }

    return gf_exp[logResult];
}

static gf* vandermonde(int nrows, int ncols) {
    int row, col, ptr;
    gf* matrix = (gf*)RS_MALLOC(nrows * ncols);
    if(NULL != matrix) {
        ptr = 0;
        for(row = 0; row < nrows; row++) {
            for(col = 0; col < ncols; col++) {
                matrix[ptr++] = galExp((gf)row, (gf)col);
            }
        }
    }

    return matrix;
}

static gf* sub_matrix(gf* matrix, int rmin, int cmin, int rmax, int cmax,  int nrows, int ncols) {
    int i, j, ptr = 0;
    gf* new_m = (gf*)RS_MALLOC( (rmax-rmin) * (cmax-cmin) );
    if(NULL != new_m) {
        for(i = rmin; i < rmax; i++) {
            for(j = cmin; j < cmax; j++) {
                new_m[ptr++] = matrix[i*ncols + j];
            }
        }
    }

    return new_m;
}

/* y = a.dot(b) */
static gf* multiply1(gf *a, int ar, int ac, gf *b, int br, int bc) {
    gf *new_m, tg;
    int r, c, i, ptr = 0;

    assert(ac == br);
    new_m = (gf*)RS_CALLOC(1, ar*bc);
    if(NULL != new_m) {

        /* this multiply is slow */
        for(r = 0; r < ar; r++) {
            for(c = 0; c < bc; c++) {
                tg = 0;
                for(i = 0; i < ac; i++) {
                    tg ^= gf_mul_table[ (a[r*ac+i] << 8) + b[i*bc+c] ];
                }

                new_m[ptr++] = tg;
            }
        }

    }

    return new_m;
}

/* copy from golang rs version */
static int code_some_shards(gf* matrixRows, gf** inputs, gf** outputs,
        int dataShards, int outputCount, int byteCount) {
    gf* in;
    int iRow, c;
    for(c = 0; c < dataShards; c++) {
        in = inputs[c];
        for(iRow = 0; iRow < outputCount; iRow++) {
            if(0 == c) {
                mul(outputs[iRow], in, matrixRows[iRow*dataShards+c], byteCount);
            } else {
                addmul(outputs[iRow], in, matrixRows[iRow*dataShards+c], byteCount);
            }
        }
    }

    return 0;
}

reed_solomon* reed_solomon_new(int data_shards, int parity_shards) {
    gf* vm = NULL;
    gf* top = NULL;
    int err = 0;
    reed_solomon* rs = NULL;

    do {
        rs = RS_MALLOC(sizeof(reed_solomon));
        if(NULL == rs) {
            return NULL;
        }
        rs->data_shards = data_shards;
        rs->parity_shards = parity_shards;
        rs->shards = (data_shards + parity_shards);
        rs->m = NULL;
        rs->parity = NULL;

        if(rs->shards > DATA_SHARDS_MAX || data_shards <= 0 || parity_shards <= 0) {
            err = 1;
            break;
        }

        vm = vandermonde(rs->shards, rs->data_shards);
        if(NULL == vm) {
            err = 2;
            break;
        }

        top = sub_matrix(vm, 0, 0, data_shards, data_shards, rs->shards, data_shards);
        if(NULL == top) {
            err = 3;
            break;
        }

        err = invert_mat(top, data_shards);
        assert(0 == err);

        rs->m = multiply1(vm, rs->shards, data_shards, top, data_shards, data_shards);
        if(NULL == rs->m) {
            err = 4;
            break;
        }

        rs->parity = sub_matrix(rs->m, data_shards, 0, rs->shards, data_shards, rs->shards, data_shards);
        if(NULL == rs->parity) {
            err = 5;
            break;
        }

        RS_FREE(vm);
        RS_FREE(top);
        vm = NULL;
        top = NULL;
        return rs;

    } while(0);

    fprintf(stderr, "err=%d\n", err);
    if(NULL != vm) {
        RS_FREE(vm);
    }
    if(NULL != top) {
        RS_FREE(top);
    }
    if(NULL != rs) {
        if(NULL != rs->m) {
            RS_FREE(rs->m);
        }
        if(NULL != rs->parity) {
            RS_FREE(rs->parity);
        }
        RS_FREE(rs);
    }

    return NULL;
}

void reed_solomon_release(reed_solomon* rs) {
    if(NULL != rs) {
        if(NULL != rs->m) {
            RS_FREE(rs->m);
        }
        if(NULL != rs->parity) {
            RS_FREE(rs->parity);
        }
        RS_FREE(rs);
    }
}

/**
 * input:
 * rs
 * data_blocks[rs->data_shards][block_size]
 * fec_blocks[rs->data_shards][block_size]
 * */
int reed_solomon_encode(reed_solomon* rs,
        unsigned char** data_blocks,
        unsigned char** fec_blocks,
        int block_size) {
    assert(NULL != rs && NULL != rs->parity);

    return code_some_shards(rs->parity, data_blocks, fec_blocks
            , rs->data_shards, rs->parity_shards, block_size);
}

/** input:
 * rs
 * original data_blocks[rs->data_shards][block_size]
 * dec_fec_blocks[nr_fec_blocks][block_size]
 * fec_block_nos: fec pos number in original fec_blocks
 * erased_blocks: erased blocks in original data_blocks
 * nr_fec_blocks: the number of erased blocks
 * */
int reed_solomon_decode(reed_solomon* rs,
        unsigned char **data_blocks,
        int block_size,
        unsigned char **dec_fec_blocks,
        unsigned int *fec_block_nos,
        unsigned int *erased_blocks,
        int nr_fec_blocks) {
    /* use stack instead of malloc, define a small number of DATA_SHARDS_MAX to save memory */
    gf dataDecodeMatrix[DATA_SHARDS_MAX*DATA_SHARDS_MAX];
    unsigned char* subShards[DATA_SHARDS_MAX];
    unsigned char* outputs[DATA_SHARDS_MAX];
    gf* m = rs->m;
    int i, j, c, swap, subMatrixRow, dataShards, nos, nshards;

    /* the erased_blocks should always sorted
     * if sorted, nr_fec_blocks times to check it
     * if not, sort it here
     * */
    for(i = 0; i < nr_fec_blocks; i++) {
        swap = 0;
        for(j = i+1; j < nr_fec_blocks; j++) {
            if(erased_blocks[i] > erased_blocks[j]) {
                /* the prefix is bigger than the following, swap */
                c = erased_blocks[i];
                erased_blocks[i] = erased_blocks[j];
                erased_blocks[j] = c;

                swap = 1;
            }
        }
        if(!swap) {
            //already sorted or sorted ok
            break;
        }
    }

    j = 0;
    subMatrixRow = 0;
    nos = 0;
    nshards = 0;
    dataShards = rs->data_shards;
    for(i = 0; i < dataShards; i++) {
        if(i != erased_blocks[j]) {
            /* this row is ok */
            for(c = 0; c < dataShards; c++) {
                dataDecodeMatrix[subMatrixRow*dataShards + c] = m[i*dataShards + c];
            }
            subShards[subMatrixRow] = data_blocks[i];
            subMatrixRow++;
        } else {
            j++;
        }
    }

    for(i = 0; i < nr_fec_blocks && subMatrixRow < dataShards; i++) {
        subShards[subMatrixRow] = dec_fec_blocks[i];
        j = dataShards + fec_block_nos[i];
        for(c = 0; c < dataShards; c++) {
            dataDecodeMatrix[subMatrixRow*dataShards + c] = m[j*dataShards + c]; //use spefic pos of original fec_blocks
        }
        subMatrixRow++;
    }

    if(subMatrixRow < dataShards) {
        //cannot correct
        return -1;
    }

    invert_mat(dataDecodeMatrix, dataShards);
    //printf("invert:\n");
    //print_matrix1(dataDecodeMatrix, dataShards, dataShards);
    //printf("nShards:\n");
    //print_matrix2(subShards, dataShards, block_size);

    for(i = 0; i < nr_fec_blocks; i++) {
        j = erased_blocks[i];
        outputs[i] = data_blocks[j];
        //data_blocks[j][0] = 0;
        memmove(dataDecodeMatrix+i*dataShards, dataDecodeMatrix+j*dataShards, dataShards);
    }
    //printf("subMatrixRow:\n");
    //print_matrix1(dataDecodeMatrix, nr_fec_blocks, dataShards);

    //printf("outputs:\n");
    //print_matrix2(outputs, nr_fec_blocks, block_size);

    return code_some_shards(dataDecodeMatrix, subShards, outputs,
            dataShards, nr_fec_blocks, block_size);
}

void test_001(void) {
    reed_solomon* rs = reed_solomon_new(11, 6);
    print_matrix1(rs->m, rs->data_shards, rs->data_shards);
    print_matrix1(rs->parity, rs->parity_shards, rs->data_shards);
}

void test_002(void) {
    char text[] = "hello world", output[256];
    int block_size = 1;
    int nrDataBlocks = sizeof(text)/sizeof(char) - 1;
    unsigned char* data_blocks[128];
    unsigned char* fec_blocks[128];
    int nrFecBlocks = 6;

    //decode
    unsigned int fec_block_nos[128], erased_blocks[128];
    unsigned char* dec_fec_blocks[128];
    int nr_fec_blocks;

    int i;
    reed_solomon* rs = reed_solomon_new(nrDataBlocks, nrFecBlocks);

    printf("%s:\n", __FUNCTION__);

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&text[i];
    }

    memset(output, 0, sizeof(output));
    memcpy(output, text, nrDataBlocks);
    for(i = 0; i < nrFecBlocks; i++) {
        fec_blocks[i] = (unsigned char*)&output[i + nrDataBlocks];
    }

    reed_solomon_encode(rs, data_blocks, fec_blocks, block_size);
    print_buf((gf*)output, "%d ", nrFecBlocks+nrDataBlocks);

    text[1] = 'x';
    text[3] = 'y';
    text[4] = 'z';
    erased_blocks[0] = 1;
    erased_blocks[1] = 3;
    erased_blocks[2] = 4;

    fec_block_nos[0] = 1;
    fec_block_nos[1] = 3;
    fec_block_nos[2] = 5;
    dec_fec_blocks[0] = fec_blocks[1];
    dec_fec_blocks[1] = fec_blocks[3];
    dec_fec_blocks[2] = fec_blocks[5];
    nr_fec_blocks = 3;

    printf("erased:%s\n", text);

    reed_solomon_decode(rs, data_blocks, block_size, dec_fec_blocks,
            fec_block_nos, erased_blocks, nr_fec_blocks);

    printf("fixed:%s\n", text);
}

int main(void) {
    fec_init();
    //test_001();
    test_002();

    return 0;
}
