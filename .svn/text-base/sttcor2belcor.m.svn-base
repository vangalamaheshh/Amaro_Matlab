function res=sttcor2belcor(prefixes,q,ignore_spins,ignore_bonds,niter,no_stt)
% average results and put in res
if ~exist('no_stt','var')
    no_stt=0;
end

for i=1:length(prefixes)
    if ~no_stt
        stt{i}=load([prefixes{i} '.stt']);
        sttlen(i)=size(stt{i},1);
    end
    cor{i}=load([prefixes{i} '.cor']);
    corlen(i)=size(cor{i},1);    
end 

if no_stt
    last=min(corlen);
else
    last=min(min(sttlen),min(corlen));
end 

if ~no_stt
    B=mean(cat(3,stt{:}),3);
    B(B>niter)=niter;
    B=B(:,(2+q*(ignore_spins)+1):end);
end

C=mean(cat(3,cor{:}),3);
res.Ts=C(:,2)';
C=C(:,(2+ignore_bonds+1):end);

for i=1:last
    if no_stt
        res.bel=[];
    else
        res.bel{i}=reshape(B(i,:)',q,size(B,2)/q);
        res.bel{i}=res.bel{i}./repmat(sum(res.bel{i}),q,1);
    end
    res.cor{i}=C(i,:)';
end 
