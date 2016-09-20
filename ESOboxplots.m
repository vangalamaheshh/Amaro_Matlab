%%%box plots for ESO barretts trios
%
%Box plots with total dels and homo dels, GD EAC vs NGD EAC
%	Episilon =.75 ; just the trios. 
ABS_seg_file=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/Aggregate_Seg_file.txt');
samples=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/TSPs_Sample_Map.txt');


ESO_samples=samples.pair_id(ismember(samples.tissue,'ESO'));
ESO_samples
ESO_seg=reorder_struct(ABS_seg_file,ismember(ABS_seg_file.sample,ESO_samples));

