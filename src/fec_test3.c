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

void test_001(void) {
    char text[] = "hello world", output[256];
    int block_size = 1;
    int nrDataBlocks = sizeof(text)/sizeof(char);
    unsigned char* data_blocks[128];
    unsigned char* fec_blocks[128];
    int nrFecBlocks = 6;

    //decode
    unsigned int fec_block_nos[128], erased_blocks[128];
    unsigned char* dec_fec_blocks[128];
    int nr_fec_blocks;

    int i;

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&text[i];
    }

    memset(output, 0, sizeof(output));
    memcpy(output, text, nrDataBlocks);
    for(i = 0; i < nrFecBlocks; i++) {
        fec_blocks[i] = (unsigned char*)&output[i + nrDataBlocks];
    }

    fec_encode(block_size, data_blocks, nrDataBlocks, fec_blocks, nrFecBlocks);
    print_buf((gf*)output, "0x%02x ", nrFecBlocks+nrDataBlocks);

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

int main(void) {
    fec_init();

    //expTable
    print_buf(gf_exp, "0x%02x ", sizeof(gf_exp)/sizeof(gf));

    //logTable
    print_int(gf_log, "%d ", sizeof(gf_log)/sizeof(int));

    //inverse
    print_buf(inverse, "0x%02x ", sizeof(inverse)/sizeof(gf));

    //mulTable
    //print_buf(gf_mul_table, "0x%02x ", sizeof(gf_mul_table)/sizeof(gf));

    printf("test001:\n");
    test_001();
    return 0;
}

