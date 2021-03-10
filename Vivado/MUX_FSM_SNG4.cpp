#include <stdbool.h>
#include "ap_int.h"

 // Include the Class
//#include <ap_shift_reg.h>

bool aa[5];
bool a[4],b[4],c[4],d[4];
bool w1[4], w2[4];
bool b1,b2,b3,b4;
ap_uint<4> temp = 1;
ap_uint<4> acc1, acc2;

int k;						//for loop control

//void MUX_SNG_3bit(ap_uint<3> *cmp, ap_uint<3> *cmp_n){
void MUX_SNG_2bit(ap_uint<2> *cmp){
	for(k=0; k<4;k++){
		aa[k] = (temp >> k) & 1;			//tokenize each bit
	}
//	aa[6] = feedback;
	aa[4] = aa[0] ^ aa[1];	//linear feedback
	temp = 0;
	for(k=0; k<4; k++){						//reconstruct register
		temp |= ap_uint<4>(aa[k+1] << k);	//shift register, record temp
		a[k] = aa[k+1];						//a = main output
	}
	if(a[3]){*cmp = 0;}
	else if(a[2]){*cmp = 1;}
	else if(a[1]){*cmp = 2;}
	else if(a[0]){*cmp = 3;}
	else{*cmp = 0;}
}
  
