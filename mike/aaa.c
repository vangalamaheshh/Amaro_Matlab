#include <string.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  if (nrhs!=1) mexErrMsgTxt("takes a single input matrix");
  if (!mxIsUint8(prhs[0])) mexErrMsgTxt("input matrix should be uint8");

  if (nlhs!=1) mexErrMsgTxt("yields a single output matrix");

  int n_inrows = mxGetM(prhs[0]);
  int n_incols = mxGetN(prhs[0]);

  printf("n_inrows %d    n_incols %d\n",n_inrows,n_incols);

  unsigned char *in = (unsigned char *)mxGetData(prhs[0]);

  int dims[2] = { n_inrows, n_incols };
  plhs[0] = mxCreateNumericArray(2,dims,mxUINT16_CLASS,mxREAL);
  unsigned short int *out = (unsigned short int *)mxGetData(plhs[0]);

  int row,col,idx,val;
  for (row=0;row<n_inrows;row++) {
    for(col=0;col<n_incols;col++) {
      idx = (col*n_inrows+row);
      val = *(in+idx);
      printf("%d\t",val);
      *(out+idx) = val+1;
    }
    printf("\n");
  }




}


