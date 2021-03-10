#include "stdio.h"
#include <string.h>
#include "math.h"
#include <stdbool.h>
#include "ap_int.h"


ap_uint<4> binary = 4;
ap_uint<2> select[30] = {2,3,0,1,2,0,0,1,0,1,\
                         0,0,0,0,1,2,3,0,1,2,\
                         0,0,1,0,1,0,0,0,0,1};
/*
ap_uint<8> binary = 60;
ap_uint<3> select[255] = {0,1,2,3,0,0,0,1,2,3,0,1,2,0,1,0,0,0,1,2,3,4,5,6,0,0,1,2,0,1,2,0,1,2,0,0,1,0,0,0,\
                          1,2,0,1,2,3,4,5,0,1,0,1,0,0,1,0,0,1,0,1,0,0,1,2,0,1,0,0,1,2,3,4,0,0,0,0,0,1,0,0,\
                          1,0,0,0,0,1,0,1,0,0,0,1,0,1,2,3,0,1,2,3,4,0,0,1,0,0,1,2,3,0,0,0,0,1,2,0,0,0,1,2,\
                          0,0,1,2,3,0,1,0,0,1,0,1,2,0,1,2,3,0,1,0,1,2,0,1,0,1,0,1,2,0,0,0,1,0,0,0,1,0,0,1,\
                          2,0,0,0,0,1,0,0,0,0,0,0,1,0,1,2,0,0,1,2,0,0,1,0,1,0,1,2,3,0,0,1,2,3,4,5,0,0,0,1,\
                          0,1,0,1,0,1,0,0,0,0,0,1,2,0,1,0,1,2,3,4,0,1,2,0,0,0,0,0,0,0,0,1,2,3,4,0,1,0,0,0,\
                          0,1,2,3,0,0,1,0,1,2,3,4,5,6,7};
*/
//unsigned char LFSR_8bit(ap_uint<8> *cmp, ap_uint<8> *cmp_n);
bool MUX_4(ap_uint<4> bin, ap_uint<2> sel, bool *sc);
//bool MUX_8(ap_uint<8> bin, ap_uint<3> sel, bool *sc);
bool output1;
int count=0;

int main(){
  
  FILE *fp;
  char *outfile = "D:\\MUX4.m";
//char *outfile = "D:\\MUX8.m";

  fp = fopen(outfile,"w");
  if (!fp) {
	  fprintf(stderr, "Can't open file %s!\r\n",outfile);
  }
  printf("File open for writing.\r\n");
  fprintf(fp,"MUX_4 SC = [\n");
//fprintf(fp,"MUX_8 SC = [\n");

  for(int i = 0; i<30; i++){
	  MUX_4(binary, select[i], &output1);
	  fprintf(fp,"%d, \n", int(output1));
	  if(int(output1)==1){
		  count++;
	  }
  }
/*
	for(int i = 0; i<255; i++){
		MUX_8(binary, select[i], &output1);
		fprintf(fp,"%d, \n", int(output1));
		if(int(output1)==1){
			count++;
		}
	}
*/
	fprintf(fp,"];");
	fclose(fp);

	printf("Sample output to file complete.\r\n");
	printf("Total 1s bit count = %d \n", count);
	return 0;
}
