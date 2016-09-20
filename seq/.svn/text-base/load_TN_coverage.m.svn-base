function C = load_TN_coverage(C)

try

C.targ = load_struct(regexprep(C.file.targ,'WITHGENES.txt$','WITHGENES_WITHGC.txt'),'%s%f%f%f%f',0);
C.targ = rename_field(C.targ,{'col1','col2','col3','col4','col5'},{'gene','chr','start','end','gc'});
C.targ = sort_struct(C.targ,{'chr','start','end'});
C.targ.len = C.targ.end-C.targ.start+1;
C.nt = slength(C.targ); C.tumcov = nan(C.nt,C.ns,C.ncat); C.normcov = nan(C.nt,C.ns,C.ncat);
for i=1:C.ns
  tmp = sort_struct(load_struct(['/xchip/tcga_scratch/lawrence/' C.sample.dir{i} '/' C.file.cov],...
    ['%s' repmat('%f',1,3+(2*C.ncat))],0),{'col2','col3','col4'});
  for c=1:C.ncat
    C.tumcov(:,i,c) = getfield(tmp,['col' num2str(4+c)]);
    C.normcov(:,i,c) = getfield(tmp,['col' num2str(4+C.ncat+c)]);
 end
end

catch me, excuse(me); end
