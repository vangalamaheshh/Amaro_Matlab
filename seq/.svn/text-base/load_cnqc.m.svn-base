function [Lt Ln R T N Z] = load_cnqc(sample,params)
% easy interface to load_and_process_CNQC_data
%
% Mike Lawrence 2009-10-29

if iscell(sample), error('Multiple samples not supported'); end
if ~exist('params','var'), params = []; end
params = impose_default_value(params,'auto',false);
params = impose_default_value(params,'flag1',2);

Lt=[];Ln=[];R=[];T=[];N=[];Z=[];

s = parse(sample,'^([^/]*)/([^/]*)/(.*)$',{'ttype','num','tech'});
if ~isempty(s), ttype=s.ttype{1}; num=s.num{1}; tech=s.tech{1};
else fprintf('  Invalid sample name.\n'); return; end
direc = ['/xchip/tcga_scratch/lawrence/' sample];
individual_name = upper(regexprep(sample,'/','-'));
tumor_lanelist = [direc '/tumor.bam.lanetable'];
normal_lanelist = [direc '/normal.bam.lanetable'];
if strcmp(ttype,'mm')
  tumor_seg = ['/xchip/tcga_scratch/lawrence/mm/analysis/20091019/mm-' num '.seg.txt'];
  normal_seg = [];
  params.log = 10;
else
  segdir = '/xchip/tcga_scratch/lawrence/db/seg';
  d = dir([segdir '/TCGA-*-' num '-0*-*-*-*.seg.data.txt']);
  if ~isempty(d), d=d(1); tumor_seg = [segdir '/' d.name]; else tumor_seg = []; end
  d = dir([segdir '/TCGA-*-' num '-1*-*-*-*.seg.data.txt']);
  if ~isempty(d), d=d(1); normal_seg = [segdir '/' d.name]; else normal_seg = []; end
end
if isempty(tumor_seg)
  d = dir([direc '/tumor*seg*txt']);
  if ~isempty(d), tumor_seg = [direc '/' d(1).name]; end
end
if isempty(normal_seg)
  d = dir([direc '/normal*seg*txt']);
  if ~isempty(d), normal_seg = [direc '/' d(1).name]; end
end

switch tech
    case 'wgs'
      tumor_rcl = [direc '/chunk1e6_cov_per_lane_tumor.txt'];
      normal_rcl = [direc '/chunk1e6_cov_per_lane_normal.txt'];
      region_list = '/xchip/tcga_scratch/lawrence/db/chunks1e6.txt';
      normals_db = '/xchip/tcga_scratch/lawrence/db/cnqc_normals/wgs/list.txt';
    case 'capture'
      if strcmp(ttype,'mm')
        tumor_rcl = [direc '/exome_cov_per_lane_tumor.txt'];
        normal_rcl = [direc '/exome_cov_per_lane_normal.txt'];
        region_list = ['/xchip/tcga_scratch/lawrence/capture/'...
          'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP_noheader.txt'];
        normals_db = '/xchip/tcga_scratch/lawrence/db/cnqc_normals/capture/list.txt';
        params.collapse_method = '10mb';
      else
        tumor_rcl = [direc '/chunk1e7_cov_per_lane_tumor.txt'];
        normal_rcl = [direc '/chunk1e7_cov_per_lane_normal.txt'];
        region_list = '/xchip/tcga_scratch/lawrence/db/chunks1e7.txt';
        switch params.flag1
          case 0, normals_db = '/xchip/tcga_scratch/lawrence/db/cnqc_normals/capture/ov_chunk1e6_list.txt';
          case 1, normals_db = '/xchip/tcga_scratch/lawrence/db/cnqc_normals/capture/ov_chunk1e7_list.txt';
          case 2, normals_db = '/xchip/tcga_scratch/lawrence/db/cnqc_normals/capture/pooled_chunk1e7_list.txt';
          otherwise, error('Unknown params.flag1');
        end
      end
    otherwise
      fprintf('  Unknown tech: should be wgs or capture\n');
      return
end

[Lt Ln R T N Z] = load_and_process_CNQC_data(tumor_rcl,normal_rcl,tumor_lanelist,normal_lanelist,...
      region_list,normals_db,tumor_seg,normal_seg,params);
