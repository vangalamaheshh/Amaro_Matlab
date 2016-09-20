function  [pon_vars,pass] = get_pon_variant_type( chr,pos,ponfile,variant_type )
% function pon_vars = get_pon_variant_type( chr,pos,ponfile,variant_type )
%   return 8 element vector of pon_info corresponding to site info vectors chr, pos
if nargin<3
    ponfile='/xchip/cga_home/petar/breast_exomeplus/pon_tokens_final/final_summed_tokens.hist.bin'
end
if nargin<4
    variant_type=repmat({'SNP'},size(chr));
end

if isempty(ponfile)
      ponfile='/xchip/cga_home/petar/breast_exomeplus/pon_tokens_final/final_summed_tokens.hist.bin'
end

% offset zero for all variant_types
p0=get_pon(chr,pos,ponfile);
pv=repmat(p0,[1,1,3]);
% DEL, INS flag
k=find(ismember(variant_type,{'DEL','INS'}));
% offset 1:2 scan for 
for o=1:2
    pv(k,:,o+1)=get_pon(chr(k),pos(k)+o,ponfile);
end
% choose 'worst-case' offset with most problems
[pbad kbad]=min(squeeze(pv(:,2,:)),[],2);

pon_vars=p0;

for i=1:3
    k1=find(kbad==i)
    pon_vars(k1,:)=pv(k1,:,i);
end

if nargout>1
    N=sum(pon_vars(1,:));
    cut1=pon_vars(:,1)>(N*0.05); %removes samples with low coverage in greater tan 5% of samples (this cutoff was pretty arbitrary, may not be necessary)
    cut2=pon_vars(:,2)<(N*0.85); %85% of PoN must have none of the other flags
    cut3=sum(pon_vars(:,3:6),2)>(N*0.05); %removes artifact events in greater than 5% of the PoN
    cut4=sum(pon_vars(:,7:8),2)>(N*0.0002); %removes germline events in greater than 0.02% of the PoN (This part was noted as a bit stringent, cuts out mutations in even one sample in PoN
    pass=~ (cut1|cut2|cut3|cut4);
end

end

function test
% $ cut -f1-17,42,66-67,82,186-187,266-268,278-279,300-311 Null_Model.maf.PoN.filtered.02-Sep-2014.maf  > Null_Model.maf.PoN.filtered.04-Sep-2014.maf
X=load_maf('/Users/stewart/Projects/Cancer/DLBCL/maf/Null_Model.filtered.maf')
% position id
X.pid=regexprep(strcat(X.Chromosome,':',cellstr(num2str(X.Start_position))),' ','');
% tabulate position id's
q=tab(X.pid)
% some sites have 30 reccurent samples 
q.x(q.n>100)
X1=trimStruct(X,ismember(X.pid,q.x(q.n>2)))
[pon_vars,pass] = get_pon_variant_type( X1.Chromosome,X1.Start_position,[],X1.Variant_Type )
XP=trimStruct(X1,pass)
tab(XP.Variant_Type)
q=tab(XP.pid)
myhist(q.n,0:20)

end