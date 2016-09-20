addpath ~/CancerGenomeAnalysis/trunk/matlab/segseq/
matdir = '/xchip/tcga_scratch2/ng/OV-0751/matfiles/';

if ~exist( 'J','var' )

    matJ = [ matdir 'OV-0751_3023J_tumor_normal_aligned_paired_reads_qual20.mat' ];
    matW = [ matdir 'OV-0751_30W1W_tumor_normal_aligned_paired_reads_qual20.mat' ];

    wigJfile = [ matdir 'OV-0751_3023J_W_400_initFP_1000_local_diff.wig.mat' ];
    wigWfile = [ matdir 'OV-0751_30W1W_W_400_initFP_1000_local_diff.wig.mat' ];

    J=load(matJ);
    W=load(matW);

    wigJ = load(wigJfile);
    wigW = load(wigWfile);
end

regionL = 102.2;
regionR = 102.27;

bkpL = 102.212495;
bkpR = 102.257878;

figure(1); clf;
subplot(2,1,1);
[nrJ,trJ]=plot_boundaries(J.READN,J.READT,4,regionL,regionR,'OV-0751 3023J');
hold on;
plot(repmat(bkpL,1,2), [0 1], 'g-','LineWidth',2);
plot(repmat(bkpR,1,2), [0 1], 'g-','LineWidth',2);
plot(repmat(102.215064,1,2), [0 1], 'r-','LineWidth',2);
plot(repmat(102.259258,1,2), [0 1], 'r-','LineWidth',2);
hold off;
subplot(2,1,2);
R4J = wigJ.R(find(wigJ.CHR==4));
POS4J = wigJ.POS(find(wigJ.CHR==4));
plot(POS4J/1e6, R4J, '.');
axis([ regionL regionR -1 1 ]);
xlabel('Chromosome 4 position (Mb)');
ylabel('Difference in log ratios');

figure(2); clf;
subplot(2,1,1);
[nrW,trW]=plot_boundaries(W.READN,W.READT,4,regionL,regionR,'OV-0751 30W1W');
hold on;
plot(repmat(bkpL,1,2), [0 1], 'g-','LineWidth',2);
plot(repmat(bkpR,1,2), [0 1], 'g-','LineWidth',2);
plot(repmat(102.214475,1,2), [0 1], 'r-','LineWidth',2);
plot(repmat(102.256744,1,2), [0 1], 'r-','LineWidth',2);
hold off;
subplot(2,1,2);
R4W = wigW.R(find(wigW.CHR==4));
POS4W = wigW.POS(find(wigW.CHR==4));
plot(POS4W/1e6, R4W, '.');
axis([ regionL regionR -1 1 ]);
xlabel('Chromosome 4 position (Mb)');
ylabel('Difference in log ratios');
