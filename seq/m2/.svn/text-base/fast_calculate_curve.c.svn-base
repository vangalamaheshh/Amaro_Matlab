/* fast_calculate_curve.c

   curve = fast_calculate_curve(mutation_positions,gene_length)
   (mutation_positions is 1-based)

   sped up version of calculate_curve.m

   To compile: (from Matlab prompt)
   >> cd /xchip/cga2/lawrence/cga/trunk/matlab/seq/m2
   >> mex fast_calculate_curve.c

*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  long n, maxpos;
  double *pos, *h;
  long i, j, ni, nj, d;

  if (nrhs!=2) mexErrMsgTxt("requires two input arguments: positions, maxpos");
  if (nlhs>1) mexErrMsgTxt("too many outputs requested");

  n = mxGetN(prhs[0]); if (n==1) n = mxGetM(prhs[0]);
  pos = mxGetPr(prhs[0]);
  maxpos = (long)mxGetScalar(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(maxpos,1,mxREAL);
  h = mxGetPr(plhs[0]);

  /* compute histogram */
  for(i=0;i<n-1;i++) {
    for(j=i+1;j<n;j++) {
      ni = (long)(*(pos+i));
      nj = (long)(*(pos+j));
      d = (ni-nj);
      if (d<0) d=-d;
      if (d<0 || d>=maxpos) continue;  /* ignore out-of-range distances */
      (*(h+d))++;
    }
  }

  /* convert to cumulative sum */
  for(i=1;i<maxpos;i++) {
    (*(h+i)) += (*(h+i-1));
  }

  /* divide by overall sum */
  for(i=0;i<maxpos;i++) {
    (*(h+i)) /= (*(h+maxpos-1));
  }

}


