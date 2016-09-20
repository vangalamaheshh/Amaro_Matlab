function aa=remove_affx_from_annot(aa)

pos=strmatch('AFFX',aa.probeset);
keep=setdiff(1:length(aa.probeset),pos);

fields={'data','probeset','refgene','desc','symb','chip'};
for i=1:length(fields)
  if isfield(aa,fields{i})
    tmp=getfield(aa,fields{i});
    if min(size(tmp))>1
      tmp=tmp(keep,:);
    else
      tmp=tmp(keep);
    end
    aa=setfield(aa,fields{i},tmp);
  end
end

%aa.data=aa.data(keep,:);
%aa.probeset=aa.probeset(keep);
%aa.refgene=aa.refgene(keep);
%aa.desc=aa.desc(keep);
%aa.symb=aa.symb(keep);

