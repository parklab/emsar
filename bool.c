#include "bool.h"

/* return value is 1 or 0 */
char getb(char* a, int index){
  int k = index/8;
  int R = 7 - index%8;
  return(a[k]>>R&1);
}

/* val is 1 or 0 */
void putb (char* a, int index, char val){
  int k = index/8;
  int R = 7 - index%8;
  a[k] = val<<R|a[k];
}

/* length of char array, given the length of the bitstring. */
int bitsize (int size){
  int newsize = size/8;
  if(size%8!=0) newsize++;
  return(newsize);
}

/* sfa_c arrays */
/* a1 : sfa_cr, s2: sfa_cp */
char getb2sfac (char* a1, char* a2, int index){
  int k= index/8;
  int R = 7 - index%8;
  if(a1[k]>>R&1) return('r');
  else if(a2[k]>>R&1) return('p');
  else return('f');
}

/* val : 'r','p','f' */
char putb2sfac (char* a1, char* a2, int index, char val){
  int k= index/8;
  int R = 7 - index%8;
  if(val=='r') { a1[k] = 1<<R|a1[k]; a2[k] = 0<<R|a2[k]; }
  else if(val=='p') { a1[k] = 0<<R|a1[k]; a2[k] = 1<<R|a2[k]; }
  else { a1[k] = 0<<R|a1[k]; a2[k] = 0<<R|a2[k]; }
}


