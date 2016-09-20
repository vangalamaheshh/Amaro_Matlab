function fh_MutSigCoverageGather(libdir,individual_set_id,mutation_list,target_list,categdir,build,build_dir,num_categories,patient_list)
% fh_MutSigCoverageGather
%
% gathers the scattered coverage jobs from ProcessCoverageForMutSig
%    1. sums the somatic_coverage.fwb files and saves a summed FWB
%    2. performs automatic category discovery
%    3. converts the somatic_coverage.txt files to MAT, and saves a summarized MAT
%
% requires the following parameters:
%    <libdir> = directory containing FixedWidthBinary.jar
%    <individual_set_id>
%    <mutation_list> = from MutSigPreprocess (so we can do automatic category discovery)
%    <target_list> = filename of list of exon targets, for coverage tabulation
%    <categdir> = directory containing all.fwb for the set of categories that was used for coverage tabulation
%    <build> = 'hg18', 'hg19', etc.
%    <build_dir> = directory containing build_info.txt, chr*.txt
%    <num_categories> = how many categories to use (set of best k categories from automatic discovery)
%                       (if this is a string beginning with "/", it will skip automatic discovery and use that file instead.)
%    <patient_list> = filename of patient list that was saved during MutSigPreprocess
%
% outputs these files:
%     <individual_set_id>.covfiles             = directory of *somatic_coverage.txt, .fwb, and .mat files
%     <individual_set_id>.coverage.mat         = compacted coverage matrix
%     <individual_set_id>.mutcategs.txt        = mutation categories (copied from categories_list)
%     <individual_set_id>.summed_coverage.fwb  = summed coverage across all samples (width-16, index=target_list)
%     <individual_set_id>.summed_coverage.fwi  = index for fwb
%
% Mike Lawrence 2010-09-28

delete_intermediate_files = true;

fprintf('fh_MutSigCoverageGather\n');
if nargin~=9
  error(['Usage: fh_MutSigCoverageGather(libdir,individual_set_id,mutation_list,'...
                      'target_list,categdir,build,build_dir,num_categories,patient_list)']);
end
demand_file(mutation_list);
demand_file(target_list);
demand_file(build_dir);
demand_file(categdir);
demand_file(patient_list);

outdir = '.';

% set java classpath to make accessible all jars in the libdir
javaclasspath([libdir;direc([libdir '/*.jar'])]);

% initialize ReferenceInfoObj from the genome_dir
try
  ReferenceInfoObj.init(build_dir);
catch me
  disp(me); disp(me.message);
  error('Failed to initialize ReferenceInfoObj from build_dir %s',build_dir);
end

% move files out of scatter job directories into covfiles directory
covdir = [outdir '/' individual_set_id '.covfiles'];
ensure_dir_exists(covdir);
system(['mv scatter*/*.somatic_coverage.* "' covdir '"']);

% gather
coverage_gather(patient_list,covdir,outdir,individual_set_id,target_list,categdir,build,build_dir,mutation_list,num_categories)

% delete intermediate files
if delete_intermediate_files
  fprintf('Deleting intermediate files\n');
  system(['rm ' covdir '/*.somatic_coverage.txt*']);
end

