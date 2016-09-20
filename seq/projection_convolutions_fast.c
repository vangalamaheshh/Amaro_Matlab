/* projection_convolutions_fast.c

   sped up version of projection_convolutions.m

   function p = projection_convolutions(Sdeg,Pdeg,score_obs,numbins,H,newH)

   input:
   Sdeg is score of each 2D degree (np,ncat+1,ncat+1)
   Pdeg is probability of each 2D degree (np,ncat+1,ncat+1)
   score_obs is observed score to calculate p-value for
   numbins is typically 1000
   H is an externally allocated matrix (numbins,1)
   newH is an externally allocated matrix (numbins,ncols)
 
   To compile: (from Matlab prompt)
   >> cd /xchip/cga2/lawrence/cga/trunk/matlab/seq
   >> mex projection_convolutions_fast.c

*/

#include <string.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
  double *Sdeg,*Pdeg;
  int *size;
  double *H, *newH;
  double score_obs;
  long numbins;
  long np,ncat,ncols,nyz,nnn;
  double binsize;
  long i,j,k,p,d1,d2,col,idx,idx2,idxa,idxb,idxc;
  double odouble;
  long o;
  double pval;
  double *returnval;

  if (nrhs!=6) mexErrMsgTxt("requires six input arguments: Sdeg,Pdeg,score_obs,numbins,H,newH");
  if (nlhs>1) mexErrMsgTxt("too many outputs requested");

  score_obs = mxGetScalar(prhs[2]);
  numbins = mxGetScalar(prhs[3]);

  size = (int *)mxGetDimensions(prhs[0]);
  np = size[0];
  ncat = size[1]-1;
  if (size[2]-1 != ncat) mexErrMsgTxt("Sdeg 2nd and 3rd dimensions should have same size");

  ncols = (ncat+1)*(ncat+2)/2;
  binsize = score_obs;
  binsize /= numbins;

  /*
  printf("np = %d\n",np);
  printf("ncat = %d\n",ncat);
  printf("ncols = %d\n",ncols);
  printf("score_obs = %g\n",score_obs);
  printf("numbins = %d\n",numbins);
  printf("binsize = %g\n",binsize);
  */

  size = (int *)mxGetDimensions(prhs[1]);
  if (size[0]!=np || size[1]-1!=ncat || size[2]-1!=ncat)
    mexErrMsgTxt("Pdeg should be same size as Sdeg");

  size = (int *)mxGetDimensions(prhs[4]);
  if (size[0]!=numbins) mexErrMsgTxt("H size should be numbins");
  
  size = (int *)mxGetDimensions(prhs[5]);
  if (size[0]!=numbins || size[1]!=ncols) mexErrMsgTxt("newH size should be (numbins,ncols)");

  Sdeg = mxGetPr(prhs[0]);
  Pdeg = mxGetPr(prhs[1]);
  H = mxGetPr(prhs[4]);
  newH = mxGetPr(prhs[5]);

  nyz = np*(ncat+1);
  nnn = numbins*ncols;

  /* initial condition: all probability is in first bin */
    for (i=1;i<numbins;i++) H[i] = 0; 
    /* memset(H,0,8*numbins); */   /* using memset is ~3x slower! */
  H[0] = 1;

  /* sequential convolution */
  for (p=0;p<np;p++) {
        for (i=0;i<nnn;i++) newH[i] = 0; 
	  /*memset(newH,0,8*nnn);*/
    col = 0;
    for (d1=0;d1<=ncat;d1++) {
      idx = p+(d1*np);
      for (d2=0;d2<=d1;d2++) {
	/*	idx = p+(d1*np)+(d2*np*(ncat+1)); */
	odouble = Sdeg[idx] / binsize;
	o = (odouble+0.5); /* round to nearest integer*/
	if (o<numbins) {
	  idx2 = (col*numbins)+o;
	  for (k=0; k<numbins-o; k++) newH[idx2+k] = Pdeg[idx] * H[k];
	} /* if */
	col++;
	idx += nyz;
      } /* next d2 */
    } /* next d1 */
    /* sum newH across columns to get H */
    for (i=0;i<numbins;i++) {
      H[i] = 0;
      idx = i;
      for (j=0;j<ncols;j++) {
	H[i] += newH[idx];
	idx += numbins;
      }
    }
  } /* next patient */

  /* calculate p-value */
  pval = 1;
  for(i=0;i<numbins;i++) pval -= H[i];

  /*
  printf("pval = %g\n",pval);
  */

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  returnval = mxGetPr(plhs[0]);
  returnval[0] = pval;


}


