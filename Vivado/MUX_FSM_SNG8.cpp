#include <stdbool.h>
#include "ap_int.h"

 // Include the Class
//#include <ap_shift_reg.h>

bool aa[9];
bool a[8],b[8];
ap_uint<8> temp = 1;
int k;						//for loop control

//void MUX_SNG_3bit(ap_uint<3> *cmp){
void MUX_SNG_3bit(ap_uint<3> *cmp, ap_uint<3> *cmp_n){
	for(k=0; k<8;k++){
		aa[k] = (temp >> k) & 1;			//tokenize each bit
	}
	aa[8] = aa[0] ^ aa[2] ^ aa[3] ^ aa[4];	//linear feedback
	temp = 0;
	for(k=0; k<8; k++){						//reconstruct register
		temp |= ap_uint<8>(aa[k+1] << k);	//shift register, record temp
		a[k] = aa[k+1];						//a = main output
		b[k] = aa[8-k];						//b = a permutated
	}

	if(a[7]){*cmp = 0;}
	else if(a[6]){*cmp = 1;}
	else if(a[5]){*cmp = 2;}
	else if(a[4]){*cmp = 3;}
	else if(a[3]){*cmp = 4;}
	else if(a[2]){*cmp = 5;}
	else if(a[1]){*cmp = 6;}
	else if(a[0]){*cmp = 7;}
	else{*cmp = 0;}

	if(b[7]){*cmp_n = 0;}
	else if(b[6]){*cmp_n = 1;}
	else if(b[5]){*cmp_n = 2;}
	else if(b[4]){*cmp_n = 3;}
	else if(b[3]){*cmp_n = 4;}
	else if(b[2]){*cmp_n = 5;}
	else if(b[1]){*cmp_n = 6;}
	else if(b[0]){*cmp_n = 7;}
	else{*cmp_n = 0;}

}
  
