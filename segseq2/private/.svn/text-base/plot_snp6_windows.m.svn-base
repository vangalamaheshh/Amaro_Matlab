function plot_snp6_windows( SNP6, RATIOS )
% FILE: plot_snp6_windows.m

%filename = '/xchip/cancergenome04/Derek/solexa/SNP6.0/HCC1143_20070801_median9_hg18.dat';
%load HCC1143_100K

filename = '/xchip/cancergenome04/Derek/solexa/SNP6.0/HCC1954_20070801_median8_hg18.dat';
%%filename = '/xchip/cancergenome04/Derek/solexa/HCC1954/HCC1954_median_hg18.dat';
%load HCC1954_100K

% Load chromosome lengths
fid = fopen('/xchip/cancergenome04/Derek/solexa/code/chromInfo_hg18.txt');
I = textscan(fid,'%u%f64%f64%s');
fclose(fid);
chrLength = I{2};
clear I;

SNP6.cn=2.^SNP6.ratios;
SNP6.cn=2.*SNP6.cn;

for c=1:23
    currChrPos = SNP6.pos( find( SNP6.chr == c ) );
    wcWindows = RATIOS.windows(find(RATIOS.chr==c));
    snp6counts = histc(currChrPos, wcWindows);
    if c==1
        windowCounts = snp6counts; 
    else
        windowCounts = [ windowCounts; snp6counts ];
    end
end
windowEnds = cumsum(windowCounts);
windowStarts = [ 1; windowEnds(1:(length(windowEnds)-1)) + 1]; 

for i=1:length(windowEnds)
    SNP6WC.chr(i) = c;
   SNP6WC.ratios(i) = median(SNP6.cn(windowStarts(i):windowEnds(i))) / 2;
   numProbes = length(windowStarts(i):windowEnds(i));
   if ( numProbes > 0 )
       SNP6WC.counts(i) = numProbes;
       SNP6WC.ci(i) = std(SNP6.cn(windowStarts(i):windowEnds(i))) / sqrt(numProbes);
   end
end

%  Plot & Summary functions
summarize=1;
if summarize==1
    selectWindows = find(isnan(SNP6WC.ratios)==0);
    plot(RATIOS.ratios(selectWindows),SNP6WC.ratios(selectWindows),'ko','MarkerSize',3,'MarkerFaceColor','k')
    hold on;
    plot(0:1:25,0:1:25,'k-')
    hold off;
    xlabel('Solexa copy ratios');
%    ylabel('238K Sty copy ratios');
    ylabel('SNP6.0 copy ratios');
    title('HCC1143 100K windows');
    
%    outfile = '/xchip/cancergenome04/Derek/solexa/SNP6.0/HCC1143_solexa_snp6.txt';
    outfile = '/xchip/cancergenome04/Derek/solexa/SNP6.0/HCC1954_solexa_snp6.txt';
%    outfile = '/xchip/cancergenome04/Derek/solexa/SNP6.0/HCC1954_solexa_sty.txt';
    fid=fopen(outfile,'w');
    % Need to FIX CI (for ratios versus log2 ratios)
    fmt=['Chromosome\tPosition\tSolexa copy\tSNP6.0 copy\r\n'];
    fprintf(fid,fmt);
    
    fmt=['%s\t%s\t%f\t%f\r\n'];
    for i=1:length(RATIOS.windows)
       fprintf(fid,fmt,num2str(RATIOS.chr(i)),num2str(RATIOS.windows(i)),...
           RATIOS.ratios(i)*2,SNP6WC.ratios(i)*2);
    end
    fclose(fid);
end
