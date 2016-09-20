function dRangerPrepare2(sample,P)
% dRangerPrepare2(sample,P)
%
% Mike Lawrence 2009-07-29

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR2');
P = impose_default_value(P,'dRangerPrepare_input_file','stringent_pairs.txt');
P = impose_default_value(P,'dRangerPrepare_output_file','stringent_pairs.mat');
P = impose_default_value(P,'mapping_quality_cutoff',5);

try

tic

for i=1:2
  if P.cancer_sample, if i==1, tn='normal';else tn='tumor'; end
  else if i==1, tn='sample';else break; end, end
  direc = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '_' P.dRangerPreprocess_output_dir_suffix];
  inputfile = [direc '/' P.dRangerPrepare_input_file];
  outputfile = [direc '/' P.dRangerPrepare_output_file];

  if ~exist(inputfile,'file'), error('%s does not exist',inputfile); end
  fprintf('Loading %s\n',inputfile);
  m = load_matrix(inputfile);
  nm = size(m,1);

  fprintf('  Processing data\n');
  % input file: (1)rgrp (2)id (3)chr1 (4)sta1 (5)end1 (6)str1 (7)qual1 (8)chr2 (9)sta2 (10)end2 (11)str2 (12)qual2 (13)flip
  % output mat: (1)id (2)chr1 (3)str1 (4)sta1 (5)end1 (6)chr2 (7)str2 (8)sta2 (9)end2 (10)qual1 (11)qual2 (12)rgrp
  x = m(:,[2 3 6 4 5 8 11 9 10 7 12 1]);

  fprintf('  Standardizing pair order\n');
  sw = (x(:,2)>x(:,6) | (x(:,2)==x(:,6) & x(:,3)>x(:,8)));
  tmp = x(sw,2:5); x(sw,2:5) = x(sw,6:9); x(sw,6:9) = tmp;
  tmp = x(sw,10); x(sw,10) = x(sw,11); x(sw,11) = tmp;
  x = [x double(sw)];  % append (13)switch

  % only filter T for mapping quality cutoff (keep all in N)  --> may want to move this step to dRangerRun
  if strcmp(tn,'tumor')
    c = P.mapping_quality_cutoff;
    x = x(x(:,10)>=c & x(:,11)>=c,:);
  end

  fprintf('  Removing duplicates\n');
  [tmp idx] = unique(x(:,2:9),'rows'); x = x(idx,:);
  nx = size(x,1);

  fprintf('  Saving MAT file (%d/%d records kept)\n',nx,nm);
  save(outputfile,'x','-v7.3');
end

% DONE
toc

catch me; excuse(me); end
