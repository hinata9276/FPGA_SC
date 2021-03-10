#include "stdio.h"
#include <string.h>
#include "math.h"
#include <stdbool.h>
#include "ap_int.h"

//unsigned char MUX_SNG_3bit(ap_uint<3> *cmp);
unsigned char MUX_SNG_3bit(ap_uint<3> *cmp, ap_uint<3> *cmp_n);
ap_uint<3> output1, output2;

int main(){
  
  FILE *fp;
  char *outfile = "F:\\8bit_MUX_FSM_SNG_output1.txt";

  fp = fopen(outfile,"w");
	if (!fp) {
		fprintf(stderr, "Can't open file %s!\r\n",outfile);
	}
	printf("File open for writing.\r\n");
	fprintf(fp,"MUX_weight_SNG = [\n");
	for(int i = 0; i<255; i++){
		MUX_SNG_3bit(&output1, &output2);
//		MUX_SNG_3bit(&output1);
//		fprintf(fp,"%d \n", int(output1));
		fprintf(fp,"%d, %d \n", int(output1), int(output2));
	}
	fprintf(fp,"];");
	fclose(fp);

	printf("Sample output to file complete.\r\n");

	return 0;
}
