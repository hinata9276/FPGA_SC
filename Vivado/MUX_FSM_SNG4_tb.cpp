#include "stdio.h"
#include <string.h>
#include "math.h"
#include <stdbool.h>
#include "ap_int.h"

//unsigned char MUX_SNG_3bit(ap_uint<3> *cmp, ap_uint<3> *cmp_n);
unsigned char MUX_SNG_2bit(ap_uint<2> *cmp);
bool ext_xor;
ap_uint<2> output1;

int main(){
  
  FILE *fp;
  char *outfile = "D:\\4bit_MUX_FSM_SNG_output.m";

  fp = fopen(outfile,"w");
	if (!fp) {
		fprintf(stderr, "Can't open file %s!\r\n",outfile);
	}
	printf("File open for writing.\r\n");
	fprintf(fp,"MUX_weight_SNG = [\n");
	for(int i = 0; i<30; i++){
		MUX_SNG_2bit(&output1);
//		MUX_SNG_3bit(&output1, &output2);
//		fprintf(fp,"%d, %d, \n", int(output1), int(output2));
		fprintf(fp,"%d,\n", int(output1));
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
