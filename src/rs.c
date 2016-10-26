#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include "fec.h"

#include "fec.c"

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

        if(rs->shards > 255 || data_shards <= 0 || parity_shards <= 0) {
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

int reed_solomon_encode(reed_solomon* rs, unsigned char** data_blocks, unsigned char** fec_blocks, int block_size) {
    return 0;
}

int reed_solomon_encode(reed_solomon* rs,
        unsigned char **data_blocks,
        int block_size, 
        unsigned char **dec_fec_blocks, 
        unsigned int *fec_block_nos,
        unsigned int *erased_blocks,
        int nr_fec_blocks) {
    return 0;
}

void test_001(void) {
    reed_solomon* rs = reed_solomon_new(11, 6);
    print_matrix1(rs->m, rs->data_shards, rs->data_shards);
    print_matrix1(rs->parity, rs->parity_shards, rs->data_shards);
}

int main(void) {
    fec_init();
    test_001();

    return 0;
}
