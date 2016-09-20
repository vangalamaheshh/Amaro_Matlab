function L = get_long_tumortype_names(S)

long = load_struct('/cga/tcga-gsc/home/lawrence/mut/analysis/long_tumor_type_names.v1.txt');
L = mapacross(S,long.short,long.long);
idx=find(strcmp('',L)); L(idx) = S(idx);
