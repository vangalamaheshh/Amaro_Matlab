function [D,v]=collapse_to_genes(C,rg,refseq_genes,rg_field_name,collapse_method)
%collapse_to_genes  Collapse SNP-granularity copy number data to gene granularity data.
%
%    [D,v]=collapse_to_genes(C,rg,refseq_genes,rg_field_name,collapse_method
%
% Collapses a D struct (with dat field m markers by n samples) to a
% new dat field (g genes by n samples).  If there are multiple
% markers covering a gene, those values are collapsed according to
% the specification in collapse_method (default 'mean') but you can
% also use 'median' or 'max','min', or 'extreme'). 
if ~exist('collapse_method','var')
  collapse_method=struct('method','mean','find_snps_type',1);
end

if isempty(rg)
  % refseq_genes is genomic locations
  m1=find(cellfun('isempty',regexp(refseq_genes,'---')));
  m2=m1;
else  
  if ~exist('rg_field_name','var') || isempty(rg_field_name)
    nms={rg(:).refseq};
  else
    nms={rg.(rg_field_name)};
  end
  if ~exist('refseq_genes','var') || isempty(refseq_genes)
      refseq_genes = unique(nms);
  end
  [M,m1,m2]=match_string_sets_hash(nms,refseq_genes);
end

D=reorder_D_rows(C,1:length(refseq_genes));
if isfield(D,'marker'), D=rmfield(D,'marker'); end;
if isfield(D,'chr'), D=rmfield(D,'chr'); end;
if isfield(D,'pos'), D=rmfield(D,'pos'); end;
if isfield(D,'cM'), D=rmfield(D,'cM'); end;
if isfield(D,'score'), D=rmfield(D,'score'); end;
if isfield(D,'cbs'), D=rmfield(D,'cbs'); end;
if isfield(D,'cbs_rl'), D=rmfield(D,'cbs_rl'); end;

D.gacc=refseq_genes;
D.gdesc=cell(length(refseq_genes),1);
D.dat=nan(length(refseq_genes),size(D.dat,2));
v=zeros(size(C.dat,1),1);
for i=1:length(m1)
  if mod(i,50)==0
    disp(i);
  end
%  chromosome2num(rg(m1(i)).chr(4:end))
  if isempty(rg)
    [chi,sti,eni]=str2genomic_location(refseq_genes{m1(i)});
    x{i}=find_snps(C,chi,sti,eni,1);
    pos_st{i}=refseq_genes{m1(i)};
  else
    x{i}=find_snps(C,chromosome2num(rg(m1(i)).chr(4:end)),rg(m1(i)).start,rg(m1(i)).end,collapse_method.find_snps_type);
    pos_st{i}=[rg(m1(i)).symb ' - ' rg(m1(i)).chr ':' num2str(rg(m1(i)).start) '-' num2str(rg(m1(i)).end)];
  end

  if ~isempty(x{i})
    v(x{i})=m2(i);
    pos_st{i}=[ pos_st{i} ' [' num2str(min(x{i})) '-' num2str(max(x{i})) ','  num2str(C.pos(min(x{i}))) '-' ...
                num2str(C.pos(max(x{i}))) ']'];
    switch collapse_method.method
     case 'mean'
      D.dat(m2(i),:)=mean(C.dat(x{i},:),1);
     case 'median'
      D.dat(m2(i),:)=median(C.dat(x{i},:),1);
     case 'min'
      D.dat(m2(i),:)=min(C.dat(x{i},:),[],1);
     case 'max'
      D.dat(m2(i),:)=max(C.dat(x{i},:),[],1);
     case 'extreme'
      min_val = min(C.dat(x{i},:),[],1);
      max_val = max(C.dat(x{i},:),[],1);
      vals = [min_val; max_val];
      [M mi] = max([abs(vals)]);
      v = zeros(1,length(min_val));
      v(find(mi==1)) = min_val(find(mi==1));
      v(find(mi==2)) = max_val(find(mi==2));
      D.dat(m2(i),:) = v;
    end
  end
end

D.gdesc(m2)=pos_st;



