function dRangerPrepare(sample,P)
% dRangerPrepare(sample,P)
%
% old style parameters: dRangerPrepare(sample,mapping_quality_cutoff)
%
% Mike Lawrence 2009

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end


P=impose_default_value(P,'cancer_sample',true);
P=impose_default_value(P,'mapping_quality_cutoff',5);

% if ~exist('mapping_quality_cutoff','var'), mapping_quality_cutoff = 5; end

try

tic

read_lengths = getBAMFileReadLength(sample,P);

id_col=1;chr1_col=2;strand1_col=3;start1_col=4;end1_col=5;chr2_col=6;strand2_col=7;
start2_col=8;end2_col=9;good1_col=10;good2_col=11;fle_col=12;switch_col=13;

for i=1:2
  if P.cancer_sample, if i==1, tn='normal';else tn='tumor'; end
  else if i==1, tn='sample';else break; end, end
  InputFile = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '_dR/all.weird.joined'];
  OutputFile = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '_dR/all.weird.joined.mat'];

  if ~exist(InputFile,'file'), error('%s does not exist',InputFile); end
  fprintf('Loading %s\n',InputFile);
  tbl = read_table(InputFile,'%f%f%f%f%f%f%f%f%f%f',char(9),0);

  fprintf('  Processing data\n');
  nw = length(tbl.dat{1}); dat = cat(2,tbl.dat{:,1:10}); clear tbl;
  % input file: (1)rgrp (2)id (3)chr1 (4)sta1 (5)str1 (6)chr2 (7)sta2 (8)str2 (9)q1 (10)q2
  % output mat: (1)id (2)chr1 (3)str1 (4)sta1 (5)end1 (6)chr2 (7)str2 (8)sta2 (9)end2 (10)q1 (11)q2 (12)rgrp
  x = [dat(:,[2 3 5 4 4 6 8 7 7 9 10 1])]; clear dat;
  x(:,[5 9]) = x(:,[5 9]) + read_lengths;
  x(:,[10 11]) = (x(:,[10 11]) >= P.mapping_quality_cutoff);

  fprintf('  Standardizing pair order\n');
  switch_flag=(x(:,chr1_col)>x(:,chr2_col) | (x(:,chr1_col)==x(:,chr2_col) & x(:,start1_col)>x(:,start2_col)));
  tmp=x(switch_flag,2:5); x(switch_flag,2:5)=x(switch_flag,6:9); x(switch_flag,6:9)=tmp;
  tmp=x(switch_flag,11); x(switch_flag,11)=x(switch_flag,10); x(switch_flag,10)=tmp;
  x=[x double(switch_flag)]; x=x(x(:,good1_col) & x(:,good2_col),:);

  fprintf('  Removing duplicates\n');
  [tmp idx] = unique(x(:,2:8),'rows'); x = x(idx,:);

  fprintf('  Saving MAT file (%d/%d records good)\n',length(x),nw);
  save(OutputFile,'x','-v7.3');
end

% DONE
toc

catch me; excuse(me); end
