function B = plot_intron_rates(G,IGR,ver,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'legend','on');
P = impose_default_value(P,'group_xlabels','size of genomic footprint');

if ~exist('ver','var')
  fprintf('Assuming ver=1');
  ver=1;
end

if ~exist('IGR','var'), IGR=[]; end

% rates by size bin
B = [];
B.min = [0   5e4 1e5 2e5 4e5 6e5 8e5 1e6]';
B.max = [5e4 1e5 2e5 4e5 6e5 8e5 1e6 inf]';

%B.min = [0   1e5 2e5 4e5 6e5 8e5]';
%B.max = [1e5 2e5 4e5 6e5 8e5 inf]';


if ~isempty(IGR), B.min(end+1)=nan; B.max(end+1)=nan; end
nbins = slength(B);
nvers = length(ver);

B.ngenes = nan(nbins,1);
B.tot_intron_N = zeros(nvers,nbins);
B.tot_intron_n = zeros(nvers,nbins);
B.rate = nan(nvers,nbins);
B.ci_low = nan(nvers,nbins);
B.ci_high = nan(nvers,nbins);

for i=1:nbins
  if ~isnan(B.min(i))
    idx = find(G.footprint>=B.min(i) & G.footprint<B.max(i));
    B.ngenes(i) = length(idx);
    for j=1:length(idx), k=idx(j);
      for v=1:nvers
        B.tot_intron_N(v,i) = B.tot_intron_N(v,i) + G.tot_intron_N{k}(ver(v));
        B.tot_intron_n(v,i) = B.tot_intron_n(v,i) + G.tot_intron_n{k}(ver(v));
      end
    end
  else  % IGR
    B.ngenes(i,1) = nan;
    for v=1:nvers
      B.tot_intron_N(v,i) = IGR.tot_N(ver(v));  B.tot_intron_n(v,i) = IGR.tot_n(ver(v));
    end
  end
  [rate ci] = binofit(B.tot_intron_n(:,i),B.tot_intron_N(:,i));
  B.rate(:,i) = 1e6*rate; B.ci_low(:,i) = 1e6*ci(:,1); B.ci_high(:,i) = 1e6*ci(:,2);
end

B.name = str2cell(sprintf('%d-%d bp (n=%d)\n',[B.min B.max B.ngenes]'));
if ~isempty(IGR), B.name{end} = 'IGR'; end
barweb(B.rate,B.ci_low,B.ci_high,0.8,P.group_xlabels);
ymax = max(B.ci_high(:)); ylim([0 1.1*ymax]);
ylabel('intron mutation rate / Mb','fontsize',20);
if strcmp(P.legend,'on'), legend(B.name,'location','eastoutside'); end
set(gca,'fontsize',20); set(gcf,'color',[1 1 1]);

return

% boxplot
G.tot_intron_rate = 1e6 * G.tot_intron_n ./ G.tot_intron_N;
x = G.tot_intron_rate(1:end-1);
g = nan(length(x),1); for i=1:length(x), g(i) = find(G.footprint(i)<B.max,1); end
figure(2); boxplot(x,g);

G2=reorder_struct(G,G.tot_intron_N/9>=0.8*G.tot_intron_len);
x = G2.tot_intron_rate(1:end-1);
g = nan(length(x),1); for i=1:length(x), g(i) = find(G2.footprint(i)<B.max,1); end
figure(3); boxplot(x,g);



% EXON rates by size bin
B = [];
B.min = [0   1e4 2e4 4e4 1e5 2e5 4e5 6e5 8e5 1e6 nan]';
B.max = [1e4 2e4 4e4 1e5 2e5 4e5 6e5 8e5 1e6 inf nan]';
nbins = slength(B);
for i=1:nbins
  if i<nbins
    idx = find(G.footprint>=B.min(i) & G.footprint<B.max(i)); B.ngenes(i,1) = length(idx);
  else  % IGR
    idx = slength(G); B.ngenes(i,1) = nan;
  end
  B.tot_exon_N(i,1) = sum(G.tot_exon_N(idx)); B.tot_exon_n(i,1) = sum(G.tot_exon_n(idx));
  [rate ci] = binofit(B.tot_exon_n(i),B.tot_exon_N(i));
  B.rate(i,1) = 1e6*rate; B.ci_low(i,1) = 1e6*ci(1); B.ci_high(i,1) = 1e6*ci(2);
end

B.name = str2cell(sprintf('%d-%d bp (n=%d)\n',[B.min B.max B.ngenes]')); B.name{end} = 'IGR';
barweb(B.rate,B.ci_low,B.ci_high,0.8,'size of genomic footprint'); ylim([0 20]);
ylabel('exon mutation rate / Mb','fontsize',20);
legend(B.name,'location','eastoutside');
set(gca,'fontsize',20); set(gcf,'color',[1 1 1]);

