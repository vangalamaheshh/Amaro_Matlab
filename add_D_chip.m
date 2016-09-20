function D=add_D_chip(D,chip_fname);

chip=read_chip_file(chip_fname);

[Mt,m1,m2]=match_string_sets_hash(chip.probeset,D.gacc);
no_match=setdiff(1:size(D.dat,1),m2);
if ~isempty(no_match)
  disp(['Did not match ' num2str(length(no_match)) ' probes']);
end

D.symb=cell(size(D.dat,1),1);
D.gdesc=cell(size(D.dat,1),1);

D.symb(m2)=chip.symb(m1);
D.gdesc(m2)=chip.desc(m1);

