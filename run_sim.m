function [R,r1,r2]=run_sim(ns,ws,delta_c,alpha,use_fdr)

prop_p=0.1;
prop_c=0.3;
prop_pc=0.1; % P(p|c)
M=10000;
M_pc=[0 M*prop_c*(1-prop_pc) M*prop_p-M*prop_c*prop_pc M*prop_c* ...
      prop_pc];
M_pc(1)=M-sum(M_pc);
sig=1;
delta_p=1;
%delta_c=1;
delta_pc=0;
%ns=[10:10:100];
%ws=[0.5:0.1:0.9];
%alpha=0.005;
%use_fdr=0;
nrep=5;
for i=1:length(ns)
  curn=ns(i);
  for j=1:length(ws)
    for n=1:nrep
      for k=1:2
        curw=ws(j);
        [D,E]=sim_data(M_pc,curn,curw,sig,delta_p,delta_c,delta_pc);
        [pv1,Br,msr,dfr,Bf,msf,dff]=ftest_many(D.dat,D.supdat(1,:),D.supdat(2,:));
        pv2=ttest2_many(D.dat,find(D.supdat(1,:)==0),find(D.supdat(1,:)==1));
        
        if (use_fdr)
          fdr1=calc_fdr_value(pv1);
          fdr2=calc_fdr_value(pv2);
        else
          fdr1=pv1;
          fdr2=pv2;
        end
        
        rej1=find(fdr1<=alpha);
        rej2=find(fdr2<=alpha);
      
        R{i,j,n,k,1}=rej1;
        R{i,j,n,k,2}=rej2; 
      end
    end        
  end        
end


for i=1:length(ns)
  curn=ns(i);
  for j=1:length(ws)
    for n=1:nrep
      for k=1:2
        curw=ws(j);

        r1(i,j,n)=length(intersect(R{i,j,n,1,1},R{i,j,n,2,1}))/ ...
                  length(union(R{i,j,n,1,1},R{i,j,n,2,1}));
        r2(i,j,n)=length(intersect(R{i,j,n,1,2},R{i,j,n,2,2}))/ ...
                  length(union(R{i,j,n,1,2},R{i,j,n,2,2}));
      end
    end
  end
end

if delta_c < 0
  dc_str=[ num2str(delta_c) 'n' ];
else
  dc_str=[ num2str(delta_c) 'p' ];
end

if use_fdr
  fdr_str='FDR';
else
  fdr_str='alpha';
end

save(['simres_' num2str(length(ns)) '_' num2str(length(ws)) '_' ...
      dc_str '_' num2str(alpha) '_' fdr_str '.mat'],'R','r1','r2');


