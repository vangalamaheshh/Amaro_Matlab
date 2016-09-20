/* 
 * Yotam Drier, yotamd@gmail.com
 */
#include "math.h"
#include "mex.h"

#define gapcost 3

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int n1, n2, m, m2, k1, k2;
    double *jumpcost;
    char *seq1;
    char *seq2;
    char *sseq;
    char *squal;
    
    signed char deltai[3] = {-1, 0, -1};
    signed char deltaj[3] = {-1, -1, 0};
    int breakpoint1=-1;
    int loconread=-1;
    int breakpoint2=-1;
    short skip=-100;
    char matchonly1=0;
    float maxsofar=0.0;
    int i, j, foriegncost;
    float max, mM1=0, mM2=0;
    unsigned long *jump, jumpcord;
    float *M1, *M2, *maxM1i;
    double *p;
    float matchscore, a, b, c, d, e;
    short endi, endj, di, dj, sofari, sofarj, *maxM1ij;
    char *which;
    
    n1 = mxGetN(prhs[0]);
    n2 = mxGetN(prhs[1]);
    m = mxGetN(prhs[2]);
    m2 = mxGetN(prhs[3]);        
    if (m!=m2) {
        printf("m=%d m2=%d\n",m,m2);
        return;
    }
    seq1 = mxCalloc(n1+1, sizeof(char));
    seq2 = mxCalloc(n2+1, sizeof(char));
    sseq = mxCalloc(m+1, sizeof(char));
    squal = mxCalloc(m+1, sizeof(char));
    mxGetString(prhs[0], seq1, n1+1);
    mxGetString(prhs[1], seq2, n2+1);
    mxGetString(prhs[2], sseq, m+1);
    mxGetString(prhs[3], squal, m+1);
    jumpcost=mxGetPr(prhs[4]);
    
    /*printf("sseq=%s squal=%s jmupcost=%f\n",sseq,squal,*jumpcost); */
    p=(double*)calloc((m+1), sizeof(double));
    M1=(float*)calloc((m+1)*(n1+1), sizeof(float));
    maxM1i=(float*)calloc(m+1, sizeof(float));
    maxM1ij=(short*)calloc(m+1, sizeof(short));
    if ((M1==0)||(maxM1i==0)||(maxM1ij==0)||(p==0)) {
        printf("memory allocation failed at 1\n");
        return;
    }
    foriegncost=gapcost+*jumpcost;
    for (i=0; i<m; i++) {
        p[i]=1-pow((double)10, -0.1*(squal[i]-33.0));
        for (j=0; j<n1; j++) {
            max=M1[i*(n1+1)+j]+(2*(seq1[j]==sseq[i])-1)*p[i];
            if (M1[(i+1)*(n1+1)+j]>M1[i*(n1+1)+j+1]) {
                if (max<M1[(i+1)*(n1+1)+j]-gapcost) {
                    max=M1[(i+1)*(n1+1)+j]-gapcost;
                }
            } else {
                if (max<M1[i*(n1+1)+j+1]-gapcost) {
                    max=M1[i*(n1+1)+j+1]-gapcost;
                }
            }
            if (max<0) max=0.0;
            M1[(i+1)*(n1+1)+j+1]=max;
            if (maxM1i[i+1]<max) {
                maxM1i[i+1]=max;
                maxM1ij[i+1]=j+1;
            }
        }
        if (mM1<maxM1i[i]) mM1=maxM1i[i];
    }
    /*printf("m=%d n1=%d mM1=%d m101=%d mij=%d %d %d %c\n%s\n",m,n1,mM1,maxM1i[100],M1[101*(n1+1)+287],M1[100*(n1+1)+286],M1[99*(n1+1)+285],seq1[286],sseq); */    
    if ((mM1<=(*jumpcost))||(mM1>=m-(*jumpcost))) {
      /*  printf("mo1 1: mM1=%f mM2=%f\n",mM1,mM2);*/
        mM2=mM1;        
        matchonly1=1;
    } else {
        /*    printf("building 2\n"); */
        M2=(float*)calloc((m+1)*(n2+1), sizeof(float));
        which=(char*)calloc((m+1)*(n2+1), sizeof(char));
        jump=(unsigned long*)calloc(m*(n2+1), sizeof(unsigned long));
        if ((M2==0)||(which==0)||(jump==0)) {
            printf("memory allocation failed at 2\n");
            return;
        }
        /*    printf("allocated\n"); */
        for (i=0; i<m; i++) {
            for (j=0; j<n2; j++) {
                matchscore=(2*(seq2[j]==sseq[i])-1)*p[i];
                a=maxM1i[i]-(*jumpcost)+matchscore;
                b=maxsofar-foriegncost+matchscore;
                c=M2[i*(n2+1)+j]+matchscore;
                d=M2[(i+1)*(n2+1)+j]-gapcost;
                e=M2[i*(n2+1)+j+1]-gapcost;
                if (a>=b) {
                    if (a>=c) {
                        if (a>=d) {
                            if (a>=e) {
                                max=a;
                                which[(i+1)*(n2+1)+j+1]=1;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        } else {
                            if (d>=e) {
                                max=d;
                                which[(i+1)*(n2+1)+j+1]=4;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        }
                    } else {
                        if (c>=d) {
                            if (c>=e) {
                                max=c;
                                which[(i+1)*(n2+1)+j+1]=3;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        } else {
                            if (d>=e) {
                                max=d;
                                which[(i+1)*(n2+1)+j+1]=4;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        }
                    }
                } else {
                    if (b>=c) {
                        if (b>=d) {
                            if (b>=e) {
                                max=b;
                                which[(i+1)*(n2+1)+j+1]=2;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        } else {
                            if (d>=e) {
                                max=d;
                                which[(i+1)*(n2+1)+j+1]=4;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        }
                    } else {
                        if (c>=d) {
                            if (c>=e) {
                                max=c;
                                which[(i+1)*(n2+1)+j+1]=3;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        } else {
                            if (d>=e) {
                                max=d;
                                which[(i+1)*(n2+1)+j+1]=4;
                            } else {
                                max=e;
                                which[(i+1)*(n2+1)+j+1]=5;
                            }
                        }
                    }
                }
                if (max<=0) {
                    max=0;
                    which[(i+1)*(n2+1)+j+1]=0;
                }
                M2[(i+1)*(n2+1)+j+1]=max;
                if (mM2<max) {
                    mM2=max;
                    endi=i+1;
                    endj=j+1;
                }
                if (which[(i+1)*(n2+1)+j+1]==1) {
                    jump[i*(n2+1)+j]=maxM1ij[i]<<10|i;
                } else {
                    if (which[(i+1)*(n2+1)+j+1]==2) {
                        jump[i*(n2+1)+j]=sofarj<<10|sofari;
                    }
                }
            }
            if (maxM1i[i]>maxsofar) {
                maxsofar=maxM1i[i];
                sofarj=maxM1ij[i];
                sofari=i;
            }
        }
        /*printf("mM2=%d endi=%d endj=%d\n",mM2,endi,endj);*/
        /*printf("done!\n");*/
        /* printf("jc=%f which=%d\n",*jumpcost,which[32*(n2+1)+156],maxM1i);
       printf("mM2=%d i=%d j=%d %d %d %d %d %d %d %d %d %d %d\n",mM2,endi,endj,M2[41*(n2+1)+165],M2[40*(n2+1)+164],M2[39*(n2+1)+163],M2[38*(n2+1)+162],M2[37*(n2+1)+161],M2[36*(n2+1)+160],M2[35*(n2+1)+159],M2[34*(n2+1)+158],M2[33*(n2+1)+157],M2[32*(n2+1)+156],M2[31*(n2+1)+155]); */
        matchonly1=mM1>=mM2;
        if (matchonly1) {
    /*        printf("mo1 2: mM1=%f mM2=%f\n",mM1,mM2);*/
            mM2=mM1;
        } else {
            for(i=endi, j=endj; (i>0)&&(j>0)&&(which[i*(n2+1)+j]>=3); i+=di, j+=dj) {
                di=deltai[which[i*(n2+1)+j]-3];
                dj=deltaj[which[i*(n2+1)+j]-3];
                /*printf("%d %d\n",di,dj);*/
            }
            /*      printf("i=%d j=%d\n",i,j);*/
            if (j==0) {
                matchonly1=1;
            } else {
                if ((i>1)&&(which[i*(n2+1)+j]>0)) {
                    jumpcord=jump[(i-1)*(n2+1)+j-1];
                    breakpoint1=jumpcord>>10;
                    breakpoint2=j-1;
                    loconread=i-1;
                    skip=i-1-(jumpcord&0x3FF);
		    if (skip==0) {		      
		      for (k1=0;(k1<loconread)&&(k1<breakpoint1)&&(seq1[breakpoint1-1-k1]!=sseq[loconread-1-k1]);k1++);
		      for (k2=0;(k2<m-loconread)&&(k2<n2-breakpoint2)&&(seq2[breakpoint2+k2]!=sseq[loconread+k2]);k2++);
		      if (k1+k2>0) {
			/*printf("k1=%d k2=%d bp1=%d bp2=%d\n",k1,k2,breakpoint1,breakpoint2);*/
			breakpoint1-=k1;
			breakpoint2+=k2;
			loconread+=k2;
			skip=(k1+k2);
			/*mM2-=foriegncost;*/
		      }
		    }
                    /* seq1(breakpoint1-loconread+1:breakpoint1)==sseq(1:loconread)
                     * //or
                     * //seq1(breakpoint1-loconread+1+skip:breakpoint1)==sseq(1:loconread-skip)
                     * //seq2(breakpoint2+1:breakpoint2+101-loconread)==sseq(loconread+1:end)    */
                } else {
                    matchonly1=1;
                }
            }
        }
        free(M2);
        free(which);
        free(jump);
    }
    free(M1);
    free(maxM1i);
    free(maxM1ij);
     
    plhs[0] = mxCreateDoubleScalar(mM2);
    plhs[1] = mxCreateLogicalScalar(matchonly1);
    plhs[2] = mxCreateDoubleScalar(breakpoint1);
    plhs[3] = mxCreateDoubleScalar(breakpoint2);
    plhs[4] = mxCreateDoubleScalar(loconread);
    plhs[5] = mxCreateDoubleScalar(skip);
    return;
}
