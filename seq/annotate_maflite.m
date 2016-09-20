function annotate_maflite(infile,outfile,build, find_DNPs)
% annotate_maflite(infile,outfile,build, find_DNPs)
%
%find_DNPs: Discover and annotate multi-nucleotides substitutions (default = true)
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
% Mike Lawrence 2009

if strcmp('GAF',build)
	fprintf('Loading GAF...');
	load /xchip/cga1/annotation/db/ucsc/hg19_GAF/GAF_F.2.mat
	annotate_maflite_GAF(infile,outfile,GAF,F, find_DNPs);
else

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

  if ~exist('build', 'var')
    error('Must provide genome build.  Can be hg17, hg18, hg18_v2, hg19, or mm9.')
  else
    P = [];
    P.build = build;
    R = load_refseq(build);
  end
  
  %%%%%%%fix to remove non-standard contigs prior to annotation

  X = reorder_struct(X,~isnan(convert_chr(X.chr)));

%  ok = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19'; ...
%  	'20';'21';'22';'X';'Y';'M';'Mt'};
%  if strncmp(build,'hg18',4)
%  	for i = 1:length(ok)
%  		ok{i} = ['chr' ok{i}];
%  	end
%  end
%  midx = ismember(X.chr, ok);
%  X = reorder_struct(X,midx);
  
  %%%%%%%%%%%%%%%%%!@ Compares both tumor alleles to ref allele to identify change
  fprintf('Preprocessing mutations: ');
  X.change = find_change_from_ref_tum1_tum2(X.start,X.end,X.ref_allele,X.tum_allele1,X.tum_allele2);
  X = reorder_struct(X, setdiff(1:slength(X),strmatch('error',X.change)));
  
  if ~exist('find_DNPs', 'var')
    find_DNPs = 1;
  else
    if ischar(find_DNPs), find_DNPs = str2double(find_DNPs); end
  end
  if find_DNPs
    fprintf('Identifying DNPS...\n');
    X = identify_DNPs_in_maf_struct(X);
  else
    fprintf('Skipping DNP identification step...\n');
  end

  
  fprintf('Classifying mutations...\n');
  X = classify_muts(X, P, R);
  fprintf('\n');
  
  fprintf('Saving output file %s\n',outfile);
  X = rmfield_if_exist(X,{'change','pos','idx'});

end
save_struct(X,outfile);

fprintf('Done!\n');
end

