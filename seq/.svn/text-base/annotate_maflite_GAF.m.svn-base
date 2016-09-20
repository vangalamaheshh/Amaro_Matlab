function annotate_maflite_GAF(infile,outfile,GAF,F,find_DNPs)
% annotate_maflite_GAF(infile,outfile,GAF,F)
%
% input format:
%
% no header line
% one line per mutation
% tab-delimited fields (all required) in this order:
%   build (e.g. '18' or 'hg18', but actually ignored)
%   chr (e.g. chr1, chr2, chrX)
%   start (1-based)
%   end (1-based)
%   ref_allele (A,C,G, or T)   or "-" for insertions
%   tum_allele1 (A,C,G, or T)  or "-" for deletions
%   tum_allele2 (A,C,G, or T)  or "-" for deletions
%   tumor_barcode
%   normal_barcode
%
% output format has header row added and the following columns appended:
%
%   gene (or "IGR" for intergenic region)
%   strand
%   classification (SNP)
%   type
%   transcript
%   genomechange
%   cDNAchange
%   codonchange
%   proteinchange
%
% GAF: UCSC GAF reference gene model struct
% F: UCSC GAF transcript sequence struct
% find_DNPs: Discover and annotate multi-nucleotides substitutions (default = true)
%
% Alex Ramos 2010

flds = {'build','chr','start','end','ref_allele','tum_allele1','tum_allele2','tumor_barcode','normal_barcode'};
%%%%%%%%%%%%%%%%%!@ Checks if file exists and format of input file
if ~exist(infile,'file'), error('Input file %s does not exist!',infile); end
d = dir(infile);
if (d.bytes==0)
  fprintf('Warning: File is empty.  Will generate a blank maf.\n');
  X = []; for i=1:length(flds), X=setfield(X,flds{i},[]); end
else
  fprintf('Loading input file %s\n',infile);
  X = load_struct(infile,'%s%s%f%f%s%s%s%s%s',char(9),0);
  if slength(X)==0, error('Input file is incorrect format'); end
  X = rename_field(X,colx(1:9),flds);
  
  if ~exist('find_DNPs', 'var')
  	find_DNPs = 1;
  end
  
  
  %%%%%%%%%%%Remove 'chr' from chromosome field if exists.
  if any(strncmp('chr', X.chr, 3))
  	for i = 1:length(X.chr)
  		if strncmp('chr', X.chr{i}, 3)
  			X.chr{i} = X.chr{i}(4:end);
  		end
  	end
  end
  
  %%%%%%%fix to remove non-standard contigs prior to annotation
  ok = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19'; ...
  	'20';'21';'22';'X';'Y';'M';'Mt'};

  midx = ismember(X.chr, ok);
  X = reorder_struct(X,midx);
  
  %%%%%%%%%%%%%%%%%!@ Compares both tumor alleles to ref allele to identify change
  fprintf('Preprocessing mutations: ');
  
  X.change = find_change_from_ref_tum1_tum2_GAF(X.start,X.end,X.ref_allele,X.tum_allele1,X.tum_allele2);
  X = reorder_struct(X, setdiff(1:slength(X),strmatch('error',X.change)));
  
  if find_DNPs
  	fprintf('Identifying DNPS...\n');
  	X = identify_DNPs_in_maf_struct(X);
  end
  
  fprintf('Classifying mutations...\n');
  X = classify_muts_GAF(X, [], GAF, F);
  fprintf('\n');
  
  fprintf('Saving output file %s\n',outfile);
  X = rmfield_if_exist(X,{'change','pos','idx'});

end
save_struct(X,outfile);

fprintf('Done!\n');

