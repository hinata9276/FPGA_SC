#include <stdbool.h>
#include "ap_int.h"

 // Include the Class
//#include <ap_shift_reg.h>

bool aa[4];
//bool aa[8];

int k;			//for loop control

//void LFSR_8bit(ap_uint<8> *cmp, ap_uint<8> *cmp_n){

void MUX_4(ap_uint<4> bin, ap_uint<2> sel, bool *sc){
	for(k=0; k<4;k++){
		aa[k] = (bin >> k) & 1;			//tokenize each bit
	}
	switch(sel){
    case 0: 
      *sc = aa[3];
      break;
    case 1:
      *sc = aa[2];
      break;
    case 2:
      *sc = aa[1];
      break;
    case 3:
      *sc = aa[0];
      break;
    default:
      break;
	}
}
/*
void MUX_8(ap_uint<8> bin, ap_uint<3> sel, bool *sc){
	for(k=0; k<8;k++){
		aa[k] = (bin >> k) & 1;			//tokenize each bit
	}
	switch(sel){
    case 0: 
      *sc = aa[7];
      break;
    case 1:
      *sc = aa[6];
      break;
    case 2:
      *sc = aa[5];
      break;
    case 3:
      *sc = aa[4];
      break;
    case 4:
      *sc = aa[3];
      break;
    case 5:
      *sc = aa[2];
      break;
    case 6:
      *sc = aa[1];
      break;
    case 7:
      *sc = aa[0];
      break;
    default:
      break;
	}
}
*/
