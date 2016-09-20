function [t n] = get_TCGA_barcodes(patient)

p = load_struct_noheader('/xchip/cga1/lawrence/ov/analysis/120patients_TCGA_barcodes.txt',...
   2,{'tumor_barcode','normal_barcode'});
p.short = regexprep(p.tumor_barcode,'TCGA-..-(\d\d\d\d)-.*','$1');
short = regexprep(patient,'.*-?(\d\d\d\d).*','$1');
t = mapacross(short,p.short,p.tumor_barcode);
n = mapacross(short,p.short,p.normal_barcode);
