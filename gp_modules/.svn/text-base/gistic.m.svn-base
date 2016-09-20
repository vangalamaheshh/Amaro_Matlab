% NOTES:
% 1) save memory
% 2) clean disks


% A new plate arrives:
% 1) Wendy runs a plate-pipeline:
%   - generate a sample info file and adds allowed normal annotation
%   - SNPFileCreator
%   - SNPtoD
%   - copy quality score and call rates based on plate normalization
%   - signal to noise metric
% 2) Wendy adds plate sample info to master 
% Q: what to do with plate updates?

% To research a set of samples:
% - select arrays for tumors and normals from master sample info file
% 2) run preprocessing 
%   - batch effect correction
%   - normal QC: (optional)
%       - normalize normals
%       - hist_qc
%       - running window median filter
%       - zero-crossing
%       - flag bad normals
%   - identify same individuals
%   - normalize
%   - hist_qc: flag bad tumors (optional)
%   - zero crossing 
%   - signal to noise metric
%   --> struct
% 3) segmentation 
%   - GLAD
%   - combine small segment
%   --> struct with segmented data
% 3.5) Remove CNVs
% 4) GISTIC core
%   - scoring
%   - permutation
%   - peel-off
%   - regions
%   --> matlab file with all the results
% 5) GISTIC output (text)
%   - write all_lesions, q-values
%   - genes_in_region
% 6) GISTIC figures


% QC of normals
%   2b) (batch correct, normalize, , hist_qc)

% Input: 
% 1) a collection of .snp files / *.mat 
% 2) a sample info file
% 3) a sample list (T+N) that you want to use
% 4) genome info file

% Module: extract_data
% combine *.snp files and extract the samples and add genomic coordinates, cyto band, and sort by location
% if needed combine platforms
% log transform

% --> mat M struct

% Module: convert M to .snp file (done)

% Module: Batch correction --> M struct

% Module: Divide_by_normals 
%   (i)   divide by mean or median of all normals (normals on plate)
%   (ii)  divide by 5 closest
%   optional divide by matched normal if exist
%   (iii) tangent
%   (iv) floor value 
% --> C struct

% Module: hist_qc 
%   (i) On normals or tumors
%   (ii) parmaters histogram
%   (iii) flag or remove
% --> C struct with flags

% Module: Filter samples based on sample info

% Module: add sample annotation

% Module: segmentation
%  (i) GLAD and parameters
%  clean temporary files 


% Module: expand segmented data and filter
%  (i) filter segmented data (number of SNPs)

% Module: remove CNVs
%  (i) table of CNVs

% Module: remove X chromosome

% Module: GISTIC statistics
%  (i) parameters, ts
% --> q p d ads score_type (in text format)

% Module: Identify regions
% --> regions

% Module: Plot GISTIC results
% --> (stats, regs)*n refseq
%  (i) ticks
%  (ii) genes
%  (iii) colors of plots
% --> pdf

% Module: Call samples
% --> regs, stats
% --> all_lesions file

% Module: List genes in regions
% --> genes, n*(sets_genes)
% --> marker for each set

%------------------------------------------------

% Module: combine platforms

% Module: Sample statistics

% Module: Gene statistics
% q-value, is it in a region, CL21 around gene

% Module: Indentify same individual

% Modules for allele-specific analysis

% Module for estimating stromal contamination

% Module to identify peaks in LOH which are not in deletion 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbstop if error
addpath ~/matlab
addpath ~/matlab/snp
addpath ~/matlab/gp_modules

cd /xchip/gistic/test

gp_gistic_preprocessing('-b','/xchip/gistic/test','-dd','/xchip/gistic/plates','-si','/xchip/gistic/test/DFCI84_070702.txt',...
                        '-o','070702', ...
                        '-of','myC.D.mat','-ns','1','-nc','1','-ncn','5','-up','1',...
                        '-bc','1','-ss','100');

% -p parameter_file -ss snp_skip -svr save_raw
