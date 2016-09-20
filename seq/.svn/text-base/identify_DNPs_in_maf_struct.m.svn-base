function F = identify_DNPs_in_maf_struct(X)

X = sort_struct(X, {'tumor_barcode','chr','start'}, [1,1,1]);  %sort in order to find DNPS

sdx = {};
sname = '';
DNPflag = 0;

%%%%%%%%finds consecutive SNPs in sorted MAF
for j=1:slength(X)
	if strncmp(X.change{j},'+',1) | strncmp(X.change{j},'-',1) | strncmp(X.change{j},'~',1)
		continue
	end

	if ~strcmp(sname, X.tumor_barcode{j}) | ~strcmp(chr, X.chr{j})
		if DNPflag, sdx = [sdx idx]; end
		sname = X.tumor_barcode{j};
		chr = X.chr{j};
		old_pos = -1;
		DNPflag = 0;
	end
	
	pos = X.start(j);
	
	if pos == (old_pos + 1)
		idx = [idx j];
		DNPflag = 1;
	elseif DNPflag == 1
		sdx = [sdx idx];
		DNPflag = 0;
		idx = [];
		idx = [idx j];
	else
		idx = [];
		idx = [idx j];
	end
	
	old_pos = pos;
end
if DNPflag, sdx = [sdx idx]; end

ff = 1:slength(X);
gg = [];
for g=1:length(sdx)
	gg = [gg sdx{g}];
end
non_DNPs = reorder_struct(X, setdiff(ff,gg)); 

%%%%%%consolidates consecutive SNP lines in MAF into single DNP lines
DNPs = {};
for k=1:length(sdx)
	Y = [];
	Yt = reorder_struct(X,sdx{k});
	Y.build = {Yt.build{1}};
	Y.chr = {Yt.chr{1}};
	Y.start = Yt.start(1);
	Y.end = Yt.end(end);

	ref = '';
	change = '';
	for h=1:slength(Yt)
		ref = [ref Yt.ref_allele{h}];
		change = [change Yt.change{h}];
	end
	
	Y.ref_allele = {ref};
	Y.tum_allele1 = {ref};
	Y.tum_allele2 = {change};
	
	Y.tumor_barcode = {Yt.tumor_barcode{1}};
	Y.normal_barcode = {Yt.normal_barcode{1}};
	
	Y.change = {['~' change]};
	
	DNPs = [DNPs Y];
end

DNPs = [DNPs non_DNPs];
F = concat_structs(DNPs);
F = sort_struct(F, {'tumor_barcode','chr','start'}, [1,1,1]); 

end