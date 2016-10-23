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

