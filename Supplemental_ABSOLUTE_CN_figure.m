seg1=load_struct('~/Downloads/TGC-Testes_DFCI_39-TP-NT-SM-7CMLJ-SM-7CMLK.segtab.txt');
seg1.Startbp=str2double(seg1.Startbp); seg1.Endbp=str2double(seg1.Endbp); seg1.modala1=str2double(seg1.modala1); seg1.modala2=str2double(seg1.modala2);
GD=0;
absolute_allelic_cn_plot( seg1,GD)

seg1=load_struct('~/Downloads/TGC-Testes_DFCI_42-TP-NT-SM-7CMLP-SM-7CMLQ.segtab.txt');
seg1.Startbp=str2double(seg1.Startbp); seg1.Endbp=str2double(seg1.Endbp); seg1.modala1=str2double(seg1.modala1); seg1.modala2=str2double(seg1.modala2);
GD=0;
absolute_allelic_cn_plot( seg1,GD)

seg1=load_struct('~/Downloads/TGC-Testes_DFCI_65-TP-NT-SM-7CMMY-SM-7CMMZ.segtab.txt');
seg1.Startbp=str2double(seg1.Startbp); seg1.Endbp=str2double(seg1.Endbp); seg1.modala1=str2double(seg1.modala1); seg1.modala2=str2double(seg1.modala2);
GD=0;
absolute_allelic_cn_plot( seg1,GD)