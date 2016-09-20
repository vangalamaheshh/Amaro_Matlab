PaperSet=load_struct('~/Documents/FullPaperSetCLL8.tsv');
Stilgenbauer=load_struct('~/Documents/AllPairsCLLStilgenbauer.tsv');
slength(Stilgenbauer)
slength(PaperSet)
StilgenbauerPaper=reorder_struct(PaperSet,ismember(PaperSet.pair_id,Stilgenbauer.pair_id));
MissingSamples=reorder_struct(Stilgenbauer,~ismember(Stilgenbauer.pair_id,StilgenbauerPaper.pair_id));

