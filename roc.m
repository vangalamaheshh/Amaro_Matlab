function [x,y,auc]=roc(pv,tr,higher_is_better)
% pv - pv (or other score). classifier is done by thresholding pv 
% if higher_is_better==1 then positives are ones with score higher
% than threshold.
% tr - vector of (0/1) indicating whether neg/pos

pv=as_column(pv);
tr=as_column(tr);

[s_pv,si]=sort(pv);
if exist('higher_is_better','var') && higher_is_better==1
  si=flipud(si);
end
s_pv=pv(si);
s_tr=tr(si);

n=length(s_tr);
tp=0;
fp=0;
fn=length(find(s_tr==1)); % all positives
tn=length(find(s_tr==0)); % all negatives


rl=runlength(s_pv');
spec=zeros(size(rl,1)+1,1);
sens=zeros(size(rl,1)+1,1);

for i=1:(size(rl,1)+1)
  if (i>1)
    np=sum(s_tr(rl(i-1,1):rl(i-1,2)));
    nn=rl(i-1,2)-rl(i-1,1)+1-np;
    tp=tp+np;
    fp=fp+nn;
    fn=fn-np;
    tn=tn-nn;
  end
  spec(i)=tn/(tn+fp+eps);
  sens(i)=tp/(tp+fn+eps);
end

x=1-spec;
y=sens;

auc=trapz(x,y);

