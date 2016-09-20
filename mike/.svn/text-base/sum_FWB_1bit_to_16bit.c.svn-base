/*    sum_FWB_1bit_to_16bit.c

      calculates sum of 1-bit data to yield 16-bit data

      given M(100000, 316)    [numbers are examples]
          where 100000 = chunk size (number of bytes from the 1-bit FWB files)
                316 = the number of FWB files

      calculates N(8*100000,1)
          where 8*100000 = number of bits from the 1-bit FWB file

      each position in N is the sum of the bits in M across the 316 samples

      Mike Lawrence 2010-12-21
*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  if (nrhs!=1) mexErrMsgTxt("takes a single input matrix");
  if (!mxIsUint8(prhs[0])) mexErrMsgTxt("input matrix should be uint8");
  if (nlhs!=1) mexErrMsgTxt("yields a single output matrix");

  int n_inrows = mxGetM(prhs[0]);
  int n_incols = mxGetN(prhs[0]);
  unsigned char *in = (unsigned char *)mxGetData(prhs[0]);

  int n_outrows = 8*n_inrows;
  int dims[2] = { n_outrows, 1 };
  plhs[0] = mxCreateNumericArray(2,dims,mxUINT16_CLASS,mxREAL);
  unsigned short int *out = (unsigned short int *)mxGetData(plhs[0]);

  unsigned char inbit = 128;
  int inrow = 0, outrow;
  for (outrow = 0; outrow < n_outrows; outrow++) {
    int inidx = inrow;
    int incol;
    for (incol = 0; incol < n_incols; incol++) {
      unsigned char inval = (unsigned char)(*(in+inidx));
      if (inval & inbit) (*(out+outrow))++;
      inidx += n_inrows;
    }
    if (inbit>1) {
      inbit >>= 1;
    } else {
      inrow++;
      inbit = 128;
    }
  }
}


