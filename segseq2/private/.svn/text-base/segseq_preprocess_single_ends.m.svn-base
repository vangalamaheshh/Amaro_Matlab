addpath /xchip/tcga_scratch/ng/segseq/code/

normalDir = '/xchip/tcga_scratch/lawrence/ov/0725/normal_SS/'
tumorDir = '/xchip/tcga_scratch/lawrence/ov/0725/tumor_SS/'

normalName = 'OV-0725_normal'
tumorName = 'OV-0725_tumor'

qualCutoff = 20

if(0)
ALLREADN = load_aligned_single_reads_dir( normalDir, qualCutoff, normalName );
ALLREADT = load_aligned_single_reads_dir( tumorDir, qualCutoff, tumorName );
end

fcName = 'OV-0725_301FF';
matfile = [ fcName '_tumor_normal_aligned_reads_qual_' num2str(qualCutoff) '.mat' ]
laneListN = 0:2
laneListT = 0:3

  N=load('/xchip/tcga_scratch/ng/OV-0725/wgs/cn/OV-0725_normal_aligned_reads_qual20.mat');
  T=load('/xchip/tcga_scratch/ng/OV-0725/wgs/cn/OV-0725_tumor_aligned_reads_qual20.mat');
ALLREADN=N.READS;
ALLREADT=T.READS;
load HG18_N36_D2_WINDOWS_100K.mat

READN = filter_reads_by_lane_num( ALLREADN, laneListN );
READT = filter_reads_by_lane_num( ALLREADT, laneListT );
RATIOS = calc_ratios_from_reads( READN, READT, WINDOWS );
save( matfile, 'READN', 'READT', 'RATIOS', '-v7.3');

clear READN;
clear READT;
clear RATIOS;

fcName = 'OV-0725_312MP';
matfile = [ fcName '_tumor_normal_aligned_reads_qual_' num2str(qualCutoff) '.mat' ]
  laneListN = 3:6
  laneListT = 4:7

load HG18_N36_D2_WINDOWS_100K.mat

READN = filter_reads_by_lane_num( ALLREADN, laneListN );
READT = filter_reads_by_lane_num( ALLREADT, laneListT );
RATIOS = calc_ratios_from_reads( READN, READT, WINDOWS );
save( matfile, 'READN', 'READT', 'RATIOS', '-v7.3');
