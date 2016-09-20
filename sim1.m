function [res,rep]=sim1(km)
if ~ispc
    addpath ~/projects/confounders
end

M_pc=[1200 600 600 300];   % M_pc=[M00 M01 M10 M11];
delta_p=0.5;
delta_c=-0.5;
delta_pc=0;
n=50;
sig=2;

use_confounder=0;

selection=struct('method','bhfdr','thresh',0.2);
test_type_conf=struct('method','ftest','confounding',[]);
test_type_noconf=struct('method','ttest');

kvec=[zeros(1,n) ones(1,0) zeros(1,n) ones(1,0)]; 


%clear res rep
for i=1:5:(n+1)
    for j=1:5:(n+1)
        if i==1
            kvec=[zeros(1,n)];
        elseif i==(n+1)
            kvec=[ones(1,n)];
        else 
            kvec=[zeros(1,n-(i-1)) ones(1,(i-1))];
        end 
        if j==1
            kvec=[kvec zeros(1,n)];
        elseif j==(n+1)
            kvec=[kvec ones(1,n)];
        else 
            kvec=[kvec zeros(1,n-(j-1)) ones(1,(j-1))];
        end 
        test_type_conf.confounding=kvec;
        for k=1:km
            [D,E]=sim_data(M_pc,n,kvec,sig,delta_p,delta_c,delta_pc);
            for c=1:2
                if c==1
                    test_type=test_type_conf;
                else
                    test_type=test_type_noconf;
                end
                [r.idx,r.q,r.p,r.s,r.pi0,r.F]=get_top_markers(D,1,test_type,-1,selection);
                        
                res{c,i,j,k}=r;
                rep(c,i,j,k,1)=1;
                rep(c,i,j,k,2)=1;
                rep(c,i,j,k,3)=1;
                
                cat=D.gsupdat(1,:)*2+D.gsupdat(2,:);
                x=res{c,i,j,k}.idx;
                if isempty(x)
                    h=zeros(1,4);
                else
                    h=histc(cat(x),0:3);
                end
                % h = m00, m01, m10, m11
                tp=h(3)+h(4);
                fp=h(1)+h(2);
                fn=M_pc(3)+M_pc(4)-tp;
                tn=M_pc(1)+M_pc(2)-fp;
                sens=tp./(tp+fn+eps);
                spec=tn./(tn+fp+eps);
                eFDR=fp./(fp+tp+eps);
                [x,y,auc]=roc(res{c,i,j,k}.F.p,[zeros(M_pc(1)+M_pc(2),1); ones(M_pc(3)+M_pc(4),1)]);
                rep(c,i,j,k,4:15)=[ h tp fp fn tn sens spec eFDR auc];
            end
        end
        disp([ i j ]);
    end
end 
