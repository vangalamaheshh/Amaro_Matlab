function [X TT NN] = fhd(individual)
fhdir = ['/xchip/cga1/firehose_output/Individual/' individual '/wgs'];
BPt = [fhdir '/ra/' individual '-Tumor.breakpoints.txt'];
BPn = [fhdir '/ra/' individual '-Normal.breakpoints.txt'];
dRmatfile =[fhdir '/ra/' individual '.dRanger_results.mat'];
[X TT NN] = load_dR_and_BP_data(individual,dRmatfile,BPt,BPn);
