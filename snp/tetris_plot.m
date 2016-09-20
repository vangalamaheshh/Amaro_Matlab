function tetris_plot(CL,tss,k,nsamples)

CL1=CL;


if k==1
  CL.dat(CL.dat<tss(1,1))=0;
else
  CL.dat(CL.dat>-tss(1,2))=0;
end

if k==1
  CL.dat=sort(CL.dat,2,'descend');
else
  CL.dat=sort(CL.dat,2,'ascend');
end

display_D(reorder_D_cols(CL,1:nsamples),[],[],'snp');
