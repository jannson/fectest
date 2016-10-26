#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include "fec.h"

#include "fec.c"

void print_buf(gf* buf, char *fmt, size_t len) {
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

void print_int(int* buf, char *fmt, size_t len) {
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

void print_matrix1(gf* matrix, int nrows, int ncols) {
    int i, j;
    for(i = 0; i < nrows; i++) {
        for(j = 0; j < ncols; j++) {
            printf("%6d ", matrix[i*ncols + j]);
        }
        printf("\n");
    }
    printf("end\n");
}

gf galExp(gf a, gf n) {
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

gf* submatrix(gf* matrix, int rmin, int cmin, int rmax, int cmax, int nrows, int ncols) {
    int i, j;
    gf* new_m = (gf*)malloc((rmax-rmin) * (cmax-cmin));
    assert(NULL != new_m);

    for(i = rmin; i < rmax; i++) {
        for(j = cmin; j < cmax; j++) {
            new_m[ (i-rmin)*(cmax-cmin) + (j-cmin) ] = matrix[i*ncols + j];
        }
    }

    return new_m;
}

gf* matrix_T(gf* b, int br, int bc) {
    gf* new_b;
    int r, c, i = 0;

    new_b = (gf*)malloc(br*bc);
    assert(NULL != new_b);

    for(c = 0; c < bc; c++) {
        for(r = 0; r < br; r++) {
            new_b[i++] = b[bc*r+c];
        }
    }

    return new_b;
}

gf* multiply1(gf* a, int ar, int ac, gf* b, int br, int bc) {
    gf *new_m, tg;
    int r, c, i;

    assert(ac == br);   //columns of a equal rows of b

    new_m = (gf*)calloc(1, ar * bc);
    assert(NULL != new_m);

    //TODO better cache
    /* new_b = matrix_T(b, br, bc);
    ptr = 0;
    for(r = 0; r < ar; r++) {
        mul(&new_m[r], &new_b[0], a[ptr++], ac);
        for(c = 1; c < ac; c++, ptr++) {
            addmul(&new_m[r], &new_b[c*br], a[ptr], ac);
        }
    } */

    for(r = 0; r < ar; r++) {
        for(c = 0; c < bc; c++) {
            tg = 0;
            for(i = 0; i < ac; i++) {
                tg ^= gf_mul_table[ (a[r*ac + i] << 8) + b[i*bc + c] ];
            }
            new_m[r*bc + c] = tg;
        }
    }

    //free(new_b);
    return new_m;
}

gf* vandermonde(int nrows, int ncols) {
    int row, col;
    gf* matrix = (gf*)calloc(1, nrows * ncols);
    assert(NULL != matrix);

    for(row = 0; row < nrows; row++) {
        for(col = 0; col < ncols; col++) {
            matrix[row*ncols + col] = galExp((gf)row, (gf)col);
        }
    }

    return matrix;
}

void codeSomeShards(gf* matrixRows, gf* intputs, gf* outputs, int dataShards, int outputCount, int byteCount) {
    gf* in;
    int r, c;

    for(c = 0; c < dataShards; c++) {
        in = &intputs[c*byteCount];
        for(r = 0; r < outputCount; r++) {
            if(0 == c) {
                mul(&outputs[r*byteCount], in, matrixRows[r*dataShards + c], byteCount);
            } else {
                addmul(&outputs[r*byteCount], in, matrixRows[r*dataShards + c], byteCount);
            }
        }
    }
}

void test_001(void) {
    printf("%s:\n", __FUNCTION__);

    //expTable
    print_buf(gf_exp, "0x%02x ", sizeof(gf_exp)/sizeof(gf));

    //logTable
    print_int(gf_log, "%d ", sizeof(gf_log)/sizeof(int));

    //inverse
    print_buf(inverse, "0x%02x ", sizeof(inverse)/sizeof(gf));

    //mulTable
    //print_buf(gf_mul_table, "0x%02x ", sizeof(gf_mul_table)/sizeof(gf));
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

    printf("%s:\n", __FUNCTION__);

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&text[i];
    }

    memset(output, 0, sizeof(output));
    memcpy(output, text, nrDataBlocks);
    for(i = 0; i < nrFecBlocks; i++) {
        fec_blocks[i] = (unsigned char*)&output[i + nrDataBlocks];
    }

    fec_encode(block_size, data_blocks, nrDataBlocks, fec_blocks, nrFecBlocks);
    print_buf((gf*)output, "%d ", nrFecBlocks+nrDataBlocks);

    text[1] = 'x';
    text[3] = 'y';
    text[4] = 'z';
    erased_blocks[0] = 1;
    erased_blocks[1] = 3;
    erased_blocks[2] = 4;

    fec_block_nos[0] = 2;
    fec_block_nos[1] = 3;
    fec_block_nos[2] = 5;
    dec_fec_blocks[0] = fec_blocks[2];
    dec_fec_blocks[1] = fec_blocks[3];
    dec_fec_blocks[2] = fec_blocks[5];
    nr_fec_blocks = 3;

    printf("erased:%s\n", text);

    fec_decode(block_size, data_blocks, nrDataBlocks,
               dec_fec_blocks, fec_block_nos, erased_blocks,
               (short)nr_fec_blocks);

    printf("fixed:%s\n", text);
}

void test_003(void) {
    int i, j;
    int a1, a2, a3, a4;

    printf("%s:\n", __FUNCTION__);

    i = 3*8;
    j = 3*5;
    a1 = gf_exp[modnn(gf_log[i] + gf_log[j]) ];

    i = 1*8;
    j = 3*5*3;
    a2 = gf_exp[modnn(gf_log[i] + gf_log[j]) ];

    i = 3*8*3;
    j = 1*5;
    a3 = gf_exp[modnn(gf_log[i] + gf_log[j]) ];

    i = 1;
    j = (3*8*3*5) % 255;
    a4 = gf_exp[modnn(gf_log[i] + gf_log[j]) ];

    printf("a1=%d a2=%d a3=%d a4=%d\n", a1, a2, a3, a4);
}

void test_004(void) {
    gf* matrix, *top, *r_m, *parity, *bshards, *outputs;
    int shards = 17;
    int dataShards = 11;
    int parityShards = 6;
    int byteCount = 1;

    //initialize
    matrix = vandermonde(shards, dataShards);
    //print_matrix1(matrix, shards, dataShards);
    //mt = matrix_T(matrix, shards, dataShards);
    //print_matrix1(mt, dataShards, shards);

    top = submatrix(matrix, 0, 0, dataShards, dataShards, shards, dataShards);
    //print_matrix1(top, dataShards, dataShards);

    invert_mat(top, dataShards);
    //print_matrix1(top, dataShards, dataShards);

    r_m = multiply1(matrix, shards, dataShards, top, dataShards, dataShards);
    printf("r_m:\n");
    print_matrix1(r_m, shards, dataShards);

    parity = submatrix(r_m, dataShards, 0, shards, dataShards, shards, dataShards);
    printf("parity:\n");
    print_matrix1(parity, parityShards, dataShards);

    //encode here
    bshards = (gf*)calloc(1, shards);
    memcpy(bshards, "hello world", dataShards); //hard code for test
    outputs = &bshards[dataShards*byteCount];
    print_buf(bshards, "%d ", dataShards);
    codeSomeShards(parity, bshards, outputs, dataShards, parityShards, byteCount);
    print_buf(bshards, "%d ", shards);

    //decode here
    {
        gf *subm, *subShards, *sub_parity;
        int i, c, shardLens[128], subMatrixRow, outputCount;

        for(i = 0; i < shards; i++) {
            shardLens[i] = byteCount;
        }

        printf("the lost is: (1:%d) (3:%d) (4:%d)\n", bshards[1], bshards[3],  bshards[4]);

        //force error
        bshards[1] = 'x';
        bshards[3] = 'y';
        bshards[4] = 'z';

        shardLens[1] = 0;
        shardLens[3] = 0;
        shardLens[4] = 0;

        subm = (gf*)malloc(dataShards * dataShards);
        subShards = (gf*)malloc(dataShards);
        subMatrixRow = 0;
        for(i = 0; i < shards && subMatrixRow < dataShards; i++) {
            if(shardLens[i] != 0) {
                for(c = 0; c < dataShards; c++) {
                    subm[subMatrixRow*dataShards + c] = r_m[i*dataShards + c];
                }

                subShards[subMatrixRow] = bshards[i];
                subMatrixRow++;
            }
        }

        invert_mat(subm, dataShards);
        printf("decodeMatrix:\n");
        print_matrix1(subm, dataShards, dataShards);
        outputs = (gf*)calloc(1, parityShards*byteCount);
        sub_parity = (gf*)calloc(1, parityShards*dataShards);
        outputCount = 0;
        for(i = 0; i < dataShards; i++) {
            if(0 == shardLens[i]) {
                //outputs[outputCount] = bshards[i]; TODO use reference
                memcpy(&sub_parity[outputCount * dataShards], &subm[i*dataShards], dataShards);
                outputCount++;
            }
        }
        printf("subMatrixRow:%d dataShards:%d outputCount:%d\n", subMatrixRow, dataShards, outputCount);
        print_matrix1(sub_parity, outputCount, dataShards);
        print_matrix1(subShards, 1, dataShards);

        //sub_parity * subShards = outputs
        codeSomeShards(sub_parity, subShards, outputs, dataShards, outputCount, byteCount); 
        printf("found lost:\n");
        print_buf(outputs, "%d ", 3);
    }
}

int main(void) {
    fec_init();

    //test_001();
    //test_002();
    //test_003();
    test_004();
    return 0;
}

