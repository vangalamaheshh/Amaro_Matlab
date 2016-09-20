function res=test_diff_anal(meta_dat,n,maxn,N,Nxv)

res=cell(2,Nxv);
cls0=1:(2*n);
cls1=(2*n)+(1:(2*n));

for i=1:Nxv
  D.dat=meta_dat((i-1)*N+(1:N),[1:(2*n) 2*maxn+(1:(2*n))]);
  verbose('1')
  [P,S]=differential_analysis(D,cls0,cls1,'ttest');
  verbose('2')
  [PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1, ...
                                                    'ttest',10000,0,10,'/xchip/data/gadgetz/lsfres/');
  verbose('3')
  [Pf,PRf]=fix_Pvalues_by_parts(P,PR,10,'/xchip/data/gadgetz/lsfres/');

  res{1,i}=Pf;
% for i=1:Nxv
%   res{1,i}=Pf((i-1)*N+(1:N));
% end
  verbose('4')
  
  [PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1, ...
                                                    'ttest',10000,1,10,'/xchip/data/gadgetz/lsfres/');
  verbose('5')
  [Pf,PRf]=fix_Pvalues_by_parts(P,PR,10,'/xchip/data/gadgetz/lsfres/');
  verbose('6')
  
  res{2,i}=Pf;
end
% for i=1:Nxv
%    res{2,i}=Pf((i-1)*N+(1:N));
% end

verbose('7')


