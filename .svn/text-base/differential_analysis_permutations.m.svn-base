function [PR,SR,OP,OS,rs]=differential_analysis_permutations(D,cls0,cls1,test_type,nperm,should_balance)
% [PR,SR,OP,OS,rs]=differential_analysis_permutations(D,cls0,cls1,test_type,nperm,should_balance)
% should_balance - should perform only balanced permutations, ones
% which approx. half are from cls0 and the remaining are from cls1. (obsolete)
%
%   CLS0 is indices in D of samples belonging to class 0.
%
%   CLS1 is indices in D of samples belonging to class 1.
%   (union(CLS0,CLs1) should equal 1:size(D.sdesc)
%
%   testtype.counfounding is a vector of length size(D.dat,2);  it contains
%   the value (bin number) of the confounder for each sample
    
if ~exist('should_balance','var')
  should_balance=1;
end

%Get number of snps, number in each class, and total number of samples
n=size(D.dat,1);
n0=length(cls0);
n1=length(cls1);
ns=n0+n1;

if should_balance
  min_n=min(n0,n1);
  if mod(min_n,2)==1
    odd_n=1;
    h0=(min_n-1)/2;
    h1=h0+1;
  else
    odd_n=0;
    h0=min_n/2;
    h1=h0;
  end
else
  cls01=[cls0 cls1];
end

% VEC is 0/1 vector giving class
vec=zeros(ns,1);
vec(cls1)=1;  

if isfield(test_type,'online')
  PR=nan(n,1);
  SR=nan(n,1);
  K=zeros(n,1);
  no_p=1;
else
  PR=nan(n,nperm);
  SR=nan(n,nperm);
  no_p=0;
end

if should_balance
  OP=PR;
  OS=SR;
else
  OP=[];
  OS=[];
end

if nperm<0
  all_perms=1;

  if should_balance
    neprm=nchoosek(n0,h0)*nchoosek(n1,h1);
    error('do not used balanced');
    %% FIX ME
  else
    if isstruct(test_type) && isfield(test_type,'confounding')
      error('all perm does not support confounding vars yet');
    else
      rs=nchoosek(1:ns,n0);
      rs=[rs zeros(size(rs,1),n1)];
      for i=1:size(rs,1)
        rs(i,n0+(1:n1))=setdiff(1:ns,rs(i,1:n0));
      end 
      nperm=size(rs,1);
    end
  end
else
  all_perms=0;
  rs=zeros(nperm,ns);
end

%% Generate Tables
% conf_all is mx1 cell array of vectors.  The ith vector contains the
% indices of the samples falling into the ith confounder category.
%
% conf_cls is mx2 table where m is number of confounding categories and
% Each (i,j) element gives the sample indices that fall into the ith
% confounder and the jth sample class.
%
% conf_all_n is an mx1 array giving the number of samples in each confounder
% class
%
% conf_cls_n is an an mx2 array giving the number of samples in each
% counfounder class and data category

if isstruct(test_type) && isfield(test_type,'confounding')
  conf_all=crosstab_sets(test_type.confounding,ones(ns,1));
  conf_cls=crosstab_sets(test_type.confounding,vec);  
  conf_all_n=cellfun('length',conf_all);
  conf_cls_n=cellfun('length',conf_cls);
end

%% 

if isfield(test_type,'gamma')
  g=test_type.gamma;
else
  g=0;
end

%% 

Binit=nperm*(n+g);
B=Binit;
active=1:n;
NP=zeros(n,1);
prm=1;
disp_step=1;
disp_delta_step=0.5; %max(50/nperm,0.005);
fprintf(1,'total permutations=%d: ',B);
diff_anal_time=0;
perm_time=0;
start_time=cputime;
loop_time=cputime;
if isfield(test_type,'report') 
    if exist(test_type.report,'file')
        delete(test_type.report);
    end
end
while ~isempty(active) && B>0
   if B<=disp_step*Binit
     fprintf(1,'%d (%d) [%d:%d,%d]...',Binit-B,length(active), ...
             cputime-start_time,diff_anal_time,perm_time);
     diff_anal_time=0;
     perm_time=0;
     disp_step=disp_step-disp_delta_step;
   end
  
  % what is the next permutation
  tmp=cputime;
  if should_balance
    error('do not balance');
    r0=randperm(n0);
    r1=randperm(n1);
    rs(prm,:)=[r0 r1];      
    perm_cls0=[ cls0(r0(1:h0)) cls1(r1(1:h1))];
    perm_cls1=[ cls0(r0(h0+(1:h1))) cls1(r1(h1+(1:h0)))];
    small_cls0=cls0(r0(1:min_n));
    small_cls1=cls1(r1(1:min_n));
    if odd_n & (rand()>0.5)
      [perm_cls0,perm_cls1]=exchange_vars(perm_cls0,perm_cls1);
    end
  else
    if all_perms
      r=rs(prm,:);
      perm_cls01=cls01(r);
      perm_cls0=perm_cls01(1:length(cls0));
      perm_cls1=perm_cls01(length(cls0)+(1:length(cls1)));
    else
      if isstruct(test_type) && isfield(test_type,'confounding')
        r=1:ns;  %vector of sample numbers
        for ci=1:size(conf_all_n,1)  %loop over all confounding bins
          rci=randperm(conf_all_n(ci));  % permute the samples in the bin
          r(conf_all{ci})=conf_all{ci}(rci); % set r so that it is a permutation of the samples, with each sample staying in its correct bin
        end
        rs(prm,:)=r; 
        perm_cls01=r(cls01);
        perm_cls0=perm_cls01(1:length(cls0));  %make the permuted classes
        perm_cls1=perm_cls01(length(cls0)+(1:length(cls1)));        
      else
        r=randperm(ns);
        rs(prm,:)=r;
        perm_cls01=cls01(r);
        perm_cls0=perm_cls01(1:length(cls0));
        perm_cls1=perm_cls01(length(cls0)+(1:length(cls1)));
      end
    end
  end
  perm_time=perm_time+cputime-tmp;

  tmp=cputime;
  if length(active)==n  %do differential analysis on permuted classes
    [cur_PR,cur_SR]=differential_analysis(D,perm_cls0, ...
                                          perm_cls1,test_type,no_p);
  else
    [cur_PR,cur_SR]=differential_analysis(reorder_D_rows(D,active),perm_cls0, ...
                                          perm_cls1,test_type,no_p);    
  end
  diff_anal_time=diff_anal_time+cputime-tmp;
  
  if length(active)==n
    if ~isempty(cur_PR)
      PR(:,prm)=cur_PR;
    end
    SR(:,prm)=cur_SR;
    NP=NP+1;
  else
    if ~isempty(cur_PR)
      PR(active,prm)=cur_PR;
    end
    SR(active,prm)=cur_SR;
    NP(active)=NP(active)+1;
  end

  if isfield(test_type,'online') && isfield(test_type,'observed')
    if length(active)==n
      ge_than_obs=find(cur_SR>=test_type.observed);
      K(ge_than_obs)=K(ge_than_obs)+1;
    else
      ge_than_obs=active(find(cur_SR>=test_type.observed(active)));
      K(ge_than_obs)=K(ge_than_obs)+1;
    end
  else
    prm=prm+1;
  end

  B=B-length(active)-g;

  if isfield(test_type,'report') 
    if isfield(test_type,'booster')
      th=test_type.booster;
    else
      th=test_type.report_th;
    end
    npa=NP(active);
    ka=K(active);    
    if isfield(test_type,'two_sided')
      ka=min(ka,npa-ka);
    end
    remove=find(2*1.96*sqrt((npa-ka+1)./((npa+3).*(ka+1)))< th);
    if ~isempty(remove)   
      f=fopen(test_type.report,'a');
      fprintf(f,'%d\t%d\t%d\t%d\t%d\t%d\n',Binit-B,length(active), ...
              length(remove),min((ka(remove)+1)./(npa(remove)+2)),cputime-loop_time,cputime-start_time);
      fclose(f);
      loop_time=cputime;
    end
  end
  
  if isfield(test_type,'booster');
    npa=NP(active);
    ka=K(active);
    if isfield(test_type,'two_sided')
      ka=min(ka,npa-ka);
    end
    keep=find(2*1.96*sqrt((npa-ka+1)./((npa+3).* ...
                                       (ka+1)))>=test_type.booster);
%    [NP(1) length(active) length(keep)]
    active=active(keep);
  end
  
  % 0 at end means return P
  if should_balance
    %    [OP(:,prm),OS(:,prm)]=differential_analysis(D,small_cls0, ...
    %                                                small_cls1,test_type);
    OP=[];
    OS=[];
  end    
end

if isfield(test_type,'online')
  PR=K;
  SR=NP;
end
fprintf(1,'\n');

%%% return as a variable
