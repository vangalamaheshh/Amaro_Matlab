function [log2ratio,PROBES] = read_dchip_file( fname, numSamples )
%  read_dchip_pos_file.m
%  INPUT:   dChipSNP file; first 3 columns contain probe info
%  OUTPUT:  cn matrix of copy numbers
%           Reads structure (R.chr, R.pos)
%
%  Derek Chiang
%  dchiang@broad.mit.edu
%


% Load coordinates of genome .marks.visual files
fid = fopen(char(fname));
C = textscan(fid,['%s%f32%f32%d8%d8' repmat('%f32',1,numSamples)],'headerLines',1);
fclose(fid);


% Assign contents to variables 
PROBES.name = C{1};
PROBES.chr = C{2};
PROBES.pos = C{3} * 1e6;

log2ratio = C{6};
for i=2:numSamples
   log2ratio = [ log2ratio C{5+i} ]; 
end
clear C;
