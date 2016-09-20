function [PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1,test_type,nperm,should_balance,nparts,lsfdir)

if nargin<7
  nparts==1;
end
if nargin<8
  lsfdir='./';
end

if nparts>1
  p=get_parts(1:nperm,nparts);
  l=lsf(lsfdir);
  h=zeros(nparts,1);
  for i=1:nparts
    [l,h(i)]=bsub(l,{'PR','SR','OP','OS','rs'},'differential_analysis_permutations',{D,cls0,cls1,test_type,length(p{i}),should_balance});
  end
  [l,res]=wait(l); % wait for all
  if isfield(test_type,'online')
    PR=zeros(size(D.dat,1),1);
    SR=zeros(size(D.dat,1),1);
    OP=PR;
    OS=SR;
    rs=zeros(nperm,length([cls0 cls1]));
    for i=1:nparts
      PR=PR+res{h(i)}.PR;
      SR=PR+res{h(i)}.SR;      
    end
    rs(p{i},:)=res{h(i)}.rs;    
  else    
    PR=zeros(size(D.dat,1),nperm);
    SR=zeros(size(D.dat,1),nperm);
    OP=PR;
    OS=SR;
    rs=zeros(nperm,length([cls0 cls1]));
    for i=1:nparts
      %    verbose(['joining results from job ' num2str(h(i)) ...
      %             ' which ran for ' num2str(res{h(i)}.end.toc)  ...
      %             ' seconds with random seed ' ...
      %             num2str(res{h(i)}.start.randseed)]);
      PR(:,p{i})=res{h(i)}.PR;
      SR(:,p{i})=res{h(i)}.SR;    
      if ~isempty(res{h(i)}.OP)
        OP(:,p{i})=res{h(i)}.OP;
      end
      if ~isempty(res{h(i)}.OS)    
        OS(:,p{i})=res{h(i)}.OS;
      end
      rs(p{i},:)=res{h(i)}.rs;    
    end
  end
else
  [PR,SR,OP,OS,rs]=differential_analysis_permutations(D,cls0,cls1,test_type,nperm,should_balance);
end
