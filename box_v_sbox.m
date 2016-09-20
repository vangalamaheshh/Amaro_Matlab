function [pmat,emat,Zmat,fmat]=box_v_sbox(box,sbox)
% sbox can include NaN fo unknown types
% box, sbox are idmat-type 
%res=((box*sbox')./max(repmat(sum(box,2),1,size(sbox,1)), ...
%		    repmat(sum(sbox,2)',size(box,1),1)));
eps=1e-6;
if nnz(isnan(sbox))>0
    nan_sbox=isnan(sbox);
    sbox(find(nan_sbox))=0;
    pmat=(box*sbox')./(repmat(sum(box,2),1,size(sbox,1))-box*nan_sbox'+eps);    
    emat=(box*sbox')./(repmat(sum(sbox,2)',size(box,1),1)+eps);
    
    C=full(sum(sbox,2))'./(size(sbox,2)-sum(nan_sbox,2)+eps)';
    B=full(sum(box,2));
    Epmat=repmat(C,size(pmat,1),1);
    spmat=sqrt(repmat(C.*(1-C),size(pmat,1),1)./(repmat(B,1,size(pmat,2))-box*nan_sbox'+eps));
    Zpmat=(pmat-Epmat)./(spmat+eps);
    %Zpmat is the same as Zemat
    Zmat=Zpmat;
else
    pmat=(box*sbox')./(repmat(sum(box,2),1,size(sbox,1))+eps);
    emat=(box*sbox')./(repmat(sum(sbox,2)',size(box,1),1)+eps);
    
    C=full(sum(sbox,2))'./(size(sbox,2)+eps);
    B=full(sum(box,2));
    Epmat=repmat(C,size(pmat,1),1);
    spmat=sqrt(repmat(C.*(1-C),size(pmat,1),1)./(repmat(B,1,size(pmat,2))+eps));
    Zpmat=(pmat-Epmat)./(spmat+eps);
    %Zpmat is the same as Zemat
    Zmat=Zpmat;
end

%X=box*sbox';
%b_sz=sum(box,2);
%sb_sz=sum(sbox,2);
fmat=[]; %zeros(size(Zmat));
%for i=1:size(fmat,1)
%    for j=1:size(fmat,2)
%       fmat(i,j)=fisher_exact_test(X(i,j),b_sz(i)-X(i,j),sb_sz(j)-X(i,j),b_sz(i)-sb_sz(j)+X(i,j));
%    end
%end
