function fh_dRanger2VCF(dRmatfile,min_somaticscoreC,ref_area,ref_build)
% fh_dRanger2VCF(individual,dRmatfile,min_somaticscore,refe_area,ref_build)
%
% Chip Stewart 2011
%
if nargin<4, error('Usage: fh_dRanger2VCF(dRmatfile,min_somaticscore,refe_area,ref_build)'); end

min_somaticscore=str2num(min_somaticscoreC);

fprintf('fh_dRanger2VCF\n');
fprintf(['  dRmatfile = ' dRmatfile '\n']);
fprintf('  min_somaticscore = %d\n', min_somaticscore);
fprintf(['  ref_area = ' ref_area '\n']);
fprintf(['  ref_build = ' ref_build '\n']);

load(dRmatfile);        

P.ref_build=ref_build;  
P.vcf_filename=regexprep(dRmatfile,'.mat$','.vcf');
P.ref_area=ref_area;
P.min_somaticscore=min_somaticscore;

V=dR2VCF4_1(X,P);
    
% skip bgzip...
% cmd=sprintf('/usr/local/bin/bgzip  %s  ',P.vcf_filename)
% system(cmd)

function test

dRmatfile='~/Projects/dRanger/data/Prostate/PR-0581.dRanger_results.detail.all.mat'
min_somaticscore=3
ref_area='/Volumes/cga1/annotation/db/ucsc/hg18_v2'
ref_build='hg18' 

fh_dRanger2VCF(dRmatfile,min_somaticscore,ref_area,ref_build)

