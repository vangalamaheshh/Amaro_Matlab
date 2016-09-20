function find_split_reads(seqs,chr1,str1,pos1,chr2,str2,pos2)

try

span1 = 2000; span2 = 2000;
min1 = pos1-round(span1/2); max1 = min1+span1-1;
min2 = pos2-round(span2/2); max2 = min2+span2-1;

ref1 = upper(genome_region(chr1,min1,max1));
if str1==1, ref1 = rc(ref1); end
ref2 = upper(genome_region(chr2,min2,max2));
if str2==0, ref2 = rc(ref2); end

% scoring matrix for extremely stringent alignment
scmat = -10*ones(size(nuc44));
for z=1:4, scmat(z,z)=1; end

seqs_all = seqs; seqs = unique(seqs); ns = length(seqs);
a = nan(ns,4); b = cell(ns,4); c = nan(2,ns,4);
for j=1:ns
  seq = seqs{j}; if strcmp(seq,'-1'), continue; end
  [a(j,1) b{j,1} c(:,j,1)] = swalign(seq,ref1,'gapopen',20,'alphabet','nt','scoringmatrix',scmat);
  [a(j,2) b{j,2} c(:,j,2)] = swalign(seq,ref2,'gapopen',20,'alphabet','nt','scoringmatrix',scmat);
  [a(j,3) b{j,3} c(:,j,3)] = swalign(rc(seq),ref1,'gapopen',20,'alphabet','nt','scoringmatrix',scmat);
  [a(j,4) b{j,4} c(:,j,4)] = swalign(rc(seq),ref2,'gapopen',20,'alphabet','nt','scoringmatrix',scmat);
  if a(j,1)+a(j,2) < a(j,3)+a(j,4)   % revcomp is better match
    seqs{j} = rc(seq);
    a(j,[1 2]) = a(j,[3 4]);
    b(j,[1 2]) = b(j,[3 4]);
    c(:,j,[1 2]) = c(:,j,[3 4]);
  end
  fprintf('(%d/%d)\t%.0f\t%.0f\n',j,ns,a(j,1),a(j,2));
end

% DIAGONAL BUBBLE PLOT
figure(1);
clf
scatter(a(:,1)+0.5*(rand(size(a,1),1)-0.5),a(:,2)+0.5*(rand(size(a,1),1)-0.5))
set(gca,'position',[0.1300    0.1100    0.7750    0.8150]);
set(gcf,'position',[1713         130         731         602]);
xlabel('length (bp) of best match to region1','fontsize',15);
ylabel('length (bp) of best match to region2','fontsize',15);
line([0 76],[76 0],'color',[0 0 0]);
set(gcf,'position',[1322          46         731         602]);

% MAKE COVERAGE PLOT
C1 = zeros(ns,span1); C2 = zeros(ns,span2);
for j=1:ns
  if strcmp(seqs{j},'-1'), continue; end
  pos = c(2,j,1);  % match to chr1
  for k=1:length(b{j,1}(3,:))
    if b{j,1}(3,k)=='-', continue; end
    C1(j,pos) = (b{j,1}(2,k)=='|'); pos=pos+1;
  end
  pos = c(2,j,2);  % match to chr2
  for k=1:length(b{j,2}(3,:))
    if b{j,2}(3,k)=='-', continue; end
    C2(j,pos) = (b{j,2}(2,k)=='|'); pos=pos+1;
  end
end

figure(2); clf
% SHOW COVERAGE PLOT
cutoff = 0; showbkpts = false;
cutoff = 15; showbkpts = false;
idx = find(a(:,1)>=cutoff & a(:,2)>=cutoff)
subplot(1,2,1); imagesc(C1(idx,:));
xlabel('position in region1','fontsize',15); ylabel('read','fontsize',15);
%xlim([900,1100]);
if showbkpts, line(1325+[0 0],[0 10000],'color',[1 1 1]); end  % actual D001 brkpt
subplot(1,2,2); imagesc(C2(idx,:));
xlabel('position in region2','fontsize',15); ylabel('read','fontsize',15);
%xlim([900,1100]);
if showbkpts, line(968+[0 0],[0 10000],'color',[1 1 1]); end  % actual D001 brkpt
set(gcf,'position',[1839         614         708         333]);

% ALIGN SPLIT READS
J = load_fasta('/xchip/tcga_scratch/lawrence/gbm/0188/wgs/D001_junction.fasta');
junc = J.seq{1}(208-80:208+80);
align_split_reads(junc,seqs);

catch me, excuse(me); end
