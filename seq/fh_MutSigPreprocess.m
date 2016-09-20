function fh_MutSigPreprocess(libdir, varargin)
% fh_MutSigPreprocess
%
% preprocesses patient list and mutations for MutSig
%
% requires the following parameters:
%    -i <individual_set_id>
%    -b <build> = 'hg18', 'hg19', etc.
%    -bd <build_dir> = directory containing <build>_info.txt and chr*.txt flatfiles for this genome build
%    -c <context_dir> = directory containing categs.txt and all.fwb to use for initial categories
%    -t <target_list> = filename of list of exon targets, for coverage tabulation
%    -cat <num_categories> = number of final categories to use (will use set of best k categories)
%                       (if this is a string beginning with "/", it will skip automatic discovery and use that file instead.)
%
% requires at least one of the following:
%    -maf1 <maffile1> = tsv file with list of maf files to load (<patient_name> <maf_filespec>)
%    -maf2 <maffile2> -maf3 <maffile3> -maf4 <maffile4> = additional maf files to combine into the final maf
%
% requires at least one of the following:
%    -wig1 <wigfile1> = tsv file with list of wig files to load (<patient_name> <wig_filespec>)
%    -wig2 <wigfile2> -wig3 <wigfile3> -wig4 <wigfile4> = additional wig files to "OR" into the final coverage calculation  
%
% optional parameters
%    -maflabel1 WGS -maflabel2 CAPTURE -maflabel3 WGS_INDEL -maflabel4 CAPTURE_INDEL = labels to write into "dataset" column
%    -p <paramfile> (optional) = file with column1=key and column2=value
%
% outputs three files:
%     <individual_set_id>.maf                     = final mutation file to use
%     <individual_set_id>.patients.txt            = final patient list to use
%     <individual_set_id>.coverage.prepare.txt    = helper file for next module (ProcessCoverageForMutSig)
%     <individual_set_id>.mutation_preprocessing_report.txt
%
% Mike Lawrence 2010-09-15

fprintf('fh_MutSigPreprocess\n');
fprintf('libdir = %s\n',libdir);
fprintf('%s\n',varargin{:});
fprintf('\n');

% set java classpath to make accessible all jars in the libdir
javaclasspath([libdir;direc([libdir '/*.jar'])]);

% FIREHOSE PARAMETERS
required_flds = {'i','c','b','bd','t','cat'};
optional_flds = {'maf1','maflabel1','maf2','maflabel2','maf3','maflabel3','maf4','maflabel4','wig1','wig2','wig3','wig4','p'};
args = handle_args([required_flds,optional_flds],varargin);
require_fields(args,required_flds);

% MUTSIG PARAMETERS
P = [];
P = process_params_file(P,args.p);

P = impose_default_value(P,'force_recalculate_maf_simplename_fields',true);
P = impose_default_value(P,'clip_tumor_from_names',true);
P = impose_default_value(P,'build',args.b);
P = impose_default_value(P,'build_dir',args.bd);
P = impose_default_value(P,'keep_barebones_fields_only',false);
P = impose_default_value(P,'consolidate_adjacent_muts_threshold',1);
P = impose_default_value(P,'remove_noncoding_mutations',true);
P = impose_default_value(P,'enforce_target_list',false);
P = impose_default_value(P,'target_list',args.t);
P = impose_default_value(P,'mutation_blacklist',[]);  % (eventually want to move these into being true FH parameters)
P = impose_default_value(P,'mutation_whitelist',[]);

disp(P)

% check categdir
categ_fwb = fullfile(args.c,'all.fwb');
categ_list = fullfile(args.c,'categs.txt');
demand_file({categ_fwb,categ_list});
tmp = load_struct(categ_list);
if strcmp(tmp.num{1},'0'), error('%s includes a zero category!',categ_file); end
mincateg = 1;
maxcateg = slength(tmp);

% initialize ReferenceInfoObj from the genome_dir
try
  ReferenceInfoObj.init(P.build_dir);
catch me
  disp(me); disp(me.message);
  error('Failed to initialize ReferenceInfoObj from P.build_dir %s',args.bd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  WIGGLE FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wigpat = [];
wigpat.name = {};
wigpat.wig = {};
for w=1:4
  wiglist = getfield(args,['wig' num2str(w)]);
  if isempty(wiglist), continue; end
  if ~exist(wiglist,'file')
    fprintf('WIG list file not found: %s\n', wiglist);
    continue;
  end
  W = load_struct_noheader(wiglist,2,{'name','wig'});
  if slength(W)==0
    fprintf('Invalid/empty WIG list file: %s\n', wiglist);
    continue;
  end
  if P.clip_tumor_from_names, W.name = regexprep(W.name,'-Tumor$',''); end
  for i=1:slength(W)
    wigfile = W.wig{i};
    if isempty(wigfile), continue; end
    if ~exist(wigfile,'file')
      fprintf('WIG file not found: %s\n', wigfile);
      continue;
    end 
    wigpat.name{end+1,1} = W.name{i};
    wigpat.wig{end+1,1} = wigfile;
  end
end
[u ui uj] = unique(wigpat.name);
wigpat2 = [];
wigpat2.name = wigpat.name(ui);
wigpat2.wig = cell(length(u),1);
for i=1:length(u), wigpat2.wig{i} = wigpat.wig(uj==i); end
wigpat = wigpat2; clear wigpat2;

if slength(wigpat)==0, error('No coverage information!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MUTATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [];
pat = [];
pat.name = {};
pat.dataset = {};
report1 = [];
for m=1:4
  maflist = getfield(args,['maf' num2str(m)]); 
  if isempty(maflist), continue; end
  if ~exist(maflist,'file')
    fprintf('MAF list file not found: %s\n', maflist);
    continue;
  end
  numread = 0;
  maflabel = getfield(args,['maflabel' num2str(m)]);
  if isempty(maflabel), maflabel = ['maf' num2str(m)]; end
  fprintf('Loading MAFs for %s\n', maflabel);
  L = load_struct_noheader(maflist,2,{'name','maf'});
  if slength(L)==0
    fprintf('Invalid/empty MAF list file: %s\n', maflist);
    continue;
  end
  if P.clip_tumor_from_names, L.name = regexprep(L.name,'-Tumor$',''); end
  for i=1:slength(L)
    if ~ismember(L.name{i},wigpat.name)
      fprintf('%s has mutation data but no coverage data: therefore is being IGNORED\n',L.name{i});
      continue;
    end
    maffile = L.maf{i};
    if isempty(maffile), continue; end
    if ~exist(maffile,'file')
      fprintf('MAF file not found: %s\n', maffile);
      continue;
    end
    try
      x{end+1,1} = load_struct(maffile);
    catch me
      fprintf('Error loading MAF file: %s\n', maffile);
      continue;
    end

    x{end} = handle_legacy_maflite_format(x{end},P); % (copies "gene" to "Hugo_Symbol", etc.)
    x{end}.firehose_patient_id = repmat({L.name{i}},slength(x{end}),1);
    x{end}.patient = x{end}.firehose_patient_id;
    x{end}.dataset = repmat({maflabel},slength(x{end}),1);

    if P.keep_barebones_fields_only
      flds = {'patient','gene','chr','start','end','type',...
        'classification','ref_allele','newbase','i_tumor_f','Protein_Change'};
      if ~isfield(x{end},'Variant_Classification')
        x{end}.Variant_Classification = repmat({''},slength(x{end}),1);
      end
      x{end} = move_to_simple_fieldnames(x{end},P);
      x{end}.newbase = find_newbase(x{end});
      x{end} = keep_fields_that_exist(x{end},flds);
      x{end}.patient = x{end}.firehose_patient_id;  % MAKE SURE
    end

    fprintf('Loaded %s\n', maffile);
    pat.name{end+1,1} = L.name{i};
    pat.dataset{end+1,1} = maflabel;
    numread = numread+1;
  end
  report1 = [report1 sprintf('Read %d MAFs of type "%s"\n',numread,maflabel)];
end

fprintf('Concatenating structs...\n');
x = concat_structs_keep_all_fields(x);

if slength(x)==0, error('No mutation data!'); end

fprintf('Calling preprocess_mutations...\n');
[x,report2] = preprocess_mutations(x,args.c,P);

if slength(x)==0, error('No mutation data survived preprocessing!'); end

fprintf('Done.\n');

% save mutation preprocessing report
report = [report1 report2];
outname = [args.i '.mutation_preprocessing_report.txt'];
save_textfile(report,outname);
fprintf('Saved %s\n', outname);

% save final maf
outname = [args.i '.maf'];
save_struct(x,outname);
fprintf('Saved %s\n', outname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PATIENT LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[u ui uj] = unique(pat.name);
pat2 = [];
pat2.name = pat.name(ui);
for i=1:length(u), pat2.dataset{i,1} = concat(pat.dataset(uj==i),'+'); end
pat=pat2; clear pat2;
fprintf('Total of %d patients had both MAF and WIG data.\n',slength(pat));
outname = [args.i '.patients.txt'];
save_struct(pat,outname);
fprintf('Saved %s\n', outname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  COVERAGE (generate "prepare" file for scatter-gather)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pat.wig = mapacross(pat.name,wigpat.name,wigpat.wig);

% scatter jobs: argument list to be passed to ExtractFromWigAllChr.java
% targetlist categfile num2str(mincateg) num2str(maxcateg) patient_name build build_dir outfile infile1 infile2 infile3 infile4

X = cell(slength(pat),1);
for i=1:slength(pat)
  X{i} = ['"' args.t '" "' categ_fwb '" ' num2str(mincateg) ' ' num2str(maxcateg) ...
          ' "' pat.name{i} '" "' P.build '" "' P.build_dir '" "' pat.name{i} '.somatic_coverage.txt"'];
  for j=1:length(pat.wig{i}), X{i} = [X{i} ' "' pat.wig{i}{j} '"']; end
end

% gather job:  argument list to be passed to fh_MutSigCoverageGather.m
% <libdir> <individual_set_id> <mutation_file> <target_list> <context_dir> <build> <build_dir> <num_categories> <patient_list_file>
% (first parameter is libdir, in order to get access to the FixedWidthBinary.jar)
% (third parameter is mutation_file, in order to perform automatic category discovery during MutSigCoverageGather)
% (last parameter is the name of the patient list file saved; its filespec will be substituted into the placeholder during prepare.py)

X{end+1} = ['"__libdir__" "' args.i '" "__mutation_list_file__" "' args.t '" "' args.c '" "' ...
         P.build '" "' P.build_dir '" "' args.cat '" "__patient_list_file__"'];

outname = [args.i '.coverage.prepare.txt'];
save_lines(X,outname);
fprintf('Saved %s\n', outname);

