%i)
SNPS=load_table('~/Documents/G201Pset5Q4.txt');
SNPS=rmfield(SNPS,'headline');
SNPS=rmfield(SNPS,'header');
plot(xhg19(SNPS.CHR,SNPS.POS,'hg18'),SNPS.OR,'b.')
ylim([0,5])
xlabel('Genome Position','FontSize',20)
ylabel('Odds Ratio','FontSize',20)

plot(xhg19(SNPS.CHR,SNPS.POS,'hg18'),-log(SNPS.P),'b.')
xlabel('Genome Position','FontSize',20)
ylabel('-log P','FontSize',20)

% ii)
% CNNM2
% These are the genes in the significantly associated region. 

%iii) His risk is 1.2272 from the odds ratio of that SNP this is also an
% extremely common variant I would say there is almost no likely hood of
% this being a causal SNP despite its strong association



