function [PR,SR,OP,OS,rs]=differential_analysis_permutations(D,cls0,cls1,test_type,nperm,should_balance)

if nargin<6
  should_balance=1;
end

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
vec=zeros(ns,1);
vec(cls1)=1;

PR=zeros(n,nperm);
SR=zeros(n,nperm);
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

if isstruct(test_type) && isfield(test_type,'confounding')
  conf_all=crosstab_sets(test_type.confounding,ones(ns,1));
  conf_cls=crosstab_sets(test_type.confounding,vec);
  conf_all_n=cellfun('length',conf_all);
  conf_cls_n=cellfun('length',conf_cls);
end

for prm=1:nperm
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
        r=1:ns;
        for ci=1:size(conf_all_n,1)
          rci=randperm(conf_all_n(ci));
          r(conf_all{ci})=conf_all{ci}(rci);
        end
        rs(prm,:)=r;
        perm_cls01=r(cls01);
        perm_cls0=perm_cls01(1:length(cls0));
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
  [PR(:,prm),SR(:,prm)]=differential_analysis(D,perm_cls0, ...
                                              perm_cls1,test_type,0); ...
      % 0 at end means return P
  if should_balance
%    [OP(:,prm),OS(:,prm)]=differential_analysis(D,small_cls0, ...
%                                                small_cls1,test_type);
    OP=[];
    OS=[];
  end    
end
