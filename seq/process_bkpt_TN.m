function [B,lanes] = process_bkpt_TN(T,N,tlanes,lnanes)

T = rename_fields(T,colx(1:8),{'fle','id','chr','start','end','strand','qual','pairmate'});
N = rename_fields(N,colx(1:8),{'fle','id','chr','start','end','strand','qual','pairmate'});
T.fle0 = T.fle - 1;   % correct to zero-based
N.fle0 = N.fle - 1;   % correct to zero-based
tlanes = load_struct('/xchip/tcga_scratch/lawrence/mm/0309/wgs/tumor_SS/lanelist.txt','%f%s',0);
nlanes = load_struct('/xchip/tcga_scratch/lawrence/mm/0309/wgs/normal_SS/lanelist.txt','%f%s',0);
tmaxlane = max(tlanes.col1); nmaxlane = max(nlanes.col1);
% combine tumor and normal
T.lane = T.fle0;
N.lane = N.fle0 + (tmaxlane+1);
tlanes.lane = tlanes.col1;
nlanes.lane = nlanes.col1 + (tmaxlane+1);
tlanes.istum = true(slength(tlanes),1);
nlanes.istum = false(slength(nlanes),1);
B = concat_structs({T,N});
lanes = concat_structs({tlanes,nlanes});
lanes = rename_fields(lanes,{'col1','col2'},{'fle','id'});
