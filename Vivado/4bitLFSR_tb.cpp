#include "stdio.h"
#include <string.h>
#include "math.h"
#include <stdbool.h>
#include "ap_int.h"

//unsigned char LFSR_8bit(ap_uint<8> *cmp, ap_uint<8> *cmp_n);
unsigned char LFSR_4bit(ap_uint<4> *cmp);
bool ext_xor;
//ap_uint<8> output1, output2;
ap_uint<4> output1;

int main(){
  
  FILE *fp;
  char *outfile = "D:\\4bitLFSRoutput.m";

  fp = fopen(outfile,"w");
	if (!fp) {
		fprintf(stderr, "Can't open file %s!\r\n",outfile);
	}
	printf("File open for writing.\r\n");
	fprintf(fp,"4 bit LFSR = [\n");
	for(int i = 0; i<30; i++){
//		LFSR_8bit(&output1, &output2);
//		fprintf(fp,"%d, %d, \n", int(output1), int(output2));
		LFSR_4bit(&output1);
		fprintf(fp,"%d, \n", int(output1));
		/*
		for(int i=0; i<8;i++){
			a[i] = (int(output1) >> i) & 1;			//tokenize each bit
		}
		ext_xor = a[0] ^ a[1];	//linear feedback
		*/
	}
	fprintf(fp,"];");
	fclose(fp);

	printf("Sample output to file complete.\r\n");

	return 0;
}
