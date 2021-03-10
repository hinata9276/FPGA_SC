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

//void LFSR_8bit(ap_uint<8> *cmp, ap_uint<8> *cmp_n){
void LFSR_4bit(ap_uint<4> *cmp){
	for(k=0; k<4;k++){
		aa[k] = (temp >> k) & 1;			//tokenize each bit
	}
//	aa[6] = feedback;
	aa[4] = aa[0] ^ aa[1];	//linear feedback
	temp = 0;
	for(k=0; k<4; k++){						//reconstruct register
		temp |= ap_uint<4>(aa[k+1] << k);	//shift register, record temp
		a[k] = aa[k+1];						//a = main output
		b[k] = !a[k];						//b = a inverted
	}
//	*cmp = temp;						//main binary output
//	*cmp_n = ~temp;						//Inverse binary output

	//////////////////////  LFSR ended here  /////////////////////////////////
	//belows are the attempt to make the 1st part of WBG work
	//so that the binary output can be shared to multiple independent SC components.
	//bool "b" array now holds the main binary output.
	//bool "c" array now holds the permutated binary output.
	//////////////////////////////////////////////////////////////////////////

	b1 = b[3] & b[2];
	w1[3] = a[3];
	w1[2] = a[2] & b[3];
	w1[1] = a[1] & b1;
	w1[0] = a[0] & b[1] & b1;

	//permutation pair
/*	for(k=0; k<8;k++){
		c[k] = a[7-k];					// c = flipped a
		d[k] = b[7-k];					// b = flipped b
	}
	b1 = d[6] & d[7];
	b2 = d[4] & d[5];
	b3 = d[3] & b2 & b1;
	b4 = d[2] & b3;
	w2[7] = c[7];						//weighted binary bit 0 for w2
	w2[6] = c[6] & d[7];				//optimized for lowest LUT (15 AND = 30 LUT/wb)
	w2[5] = c[5] & b1;
	w2[4] = c[4] & d[5] & b1;
	w2[3] = c[3] & b2 & b1;
	w2[2] = c[2] & b3;
	w2[1] = c[1] & b4;
	w2[0] = c[0] & d[1] & b4;
*/
	acc1 = 0;
//	acc2 = 0;							//permutated pair
	for(k=0; k<4; k++){
		acc1 |= ap_uint<4>(w1[k] << k);		//concatenate w1 bits into byte
//		acc2 |= ap_uint<8>(w2[k] << k);		//concatenate w2 bits into byte (permutated pair)
		}
	*cmp = acc1;						//wb output 1, to be shared among other wbgs
//	*cmp_n = acc2;						//wb output 2, to be shared among other wbgs (permutated pair)
}
  
