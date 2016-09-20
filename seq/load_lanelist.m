function L = load_lanelist(sample,tn)

if lower(tn(1))=='t', tn='tumor';
elseif lower(tn(1))=='n', tn='normal';
else error('"tn" should be "tumor" or "normal"'); end

file = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '.bam.lanelist'];
L = rename_fields(load_struct(file,'%f%s',0),{'col1','col2'},{'lane','id'});
