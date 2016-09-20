function [X,R] = analyze_454_rearrangement_reads(X,R,P)
% [X,R] = analyze_454_rearrangement_reads(X,R,P)
%
% X is struct of 454 reads from sffread + convert_sff
% R is dRanger validation design file
%
% Mike Lawrence 2010-04-20

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'outdir','*required*');
P = impose_default_value(P,'range_to_process',[]);
P = impose_default_value(P,'min_readlength',200);
P = impose_default_value(P,'leader_sequences','TCAGAGACAC|TCAGAGTTGGC');
 these sequences need to be determined for each region

require_fields(X,{'region','header','seq','qual','len','trim'});
require_fields(R,{'amplicon','leftprimer','rightprimer'}); %,'chr1','str1','pos1','chr2','str2','pos2'});

if ~exist(P.outdir,'dir'), mkdir(P.outdir); end

X = reorder_struct(X,X.trim>=P.min_readlength);
nx = slength(X);
nr = slength(R);

for i=1:nx, X.seq{i} = upper(X.seq{i}(1:X.trim(i))); X.qual{i}=X.qual{i}(1:X.trim(i)); end

X.seq = regexprep(X.seq,['^(' P.leader_sequences ')'],'');

% readlength histograms
[u ui uj] = unique(X.region);
for r=1:length(u)
  hist(X.trim(X.region==u(r)));
  xlim([min(X.trim) max(X.trim)]);
  print_to_file([P.outdir '/readlength_histogram_region' num2str(u(r)) '.png']);
end

R.leftprimer = upper(R.leftprimer);
R.rightprimer = upper(R.rightprimer);
R.amplicon = upper(R.amplicon);
%R.amplicon = regexprep(upper(R.amplicon),'N+','NNN');

nlp = cellfun('length',R.leftprimer);
nrp = cellfun('length',R.rightprimer);
mnp = min([nlp;nrp]);
if min(X.trim)<mnp, error('some read lengths are shorter than shortest primer'); end

% for each read, see if 5' end (after removing leader sequence)
%    is an exact match to first N=mnp bases of any leftprimer or rightprimer

matchL = false(nx,nr);
matchR = false(nx,nr);
fprintf('Assigning to amplicons: ');
for x=1:nx, if ~mod(x,10000), fprintf('%d/%d ',x,nx); end
  matchL(x,:) = strncmp(X.seq{x}(1:mnp),R.leftprimer,mnp);
  matchR(x,:) = strncmp(X.seq{x}(1:mnp),R.rightprimer,mnp);
end,fprintf('\n');

% for each read, find out which amplicon it comes from (if unambiguously determined)

X.amp = zeros(nx,1);
for x=1:nx
  idx = find(matchL(x,:)|matchR(x,:));
  if length(idx)==1, X.amp(x)=idx; end
end

% amplicon hit histograms
figure(1);clf;hold on
[u ui uj] = unique(X.region);
for r=1:length(u)
  subplot(length(u),1,r);
  hist(X.amp(X.region==u(r)),1:nr),xlim([1 nr]),ylim([0 1000]),line([192.5 192.5],[0 1000]);
end
hold off
print_to_file([P.outdir '/hits.png']);

% count alignments to each amplicon
X.hit = zeros(nx,1);

scmat = -2*ones(size(nuc44)); for z=1:4, scmat(z,z)=4; end
swparams = {'gapopen',32,'extendgap',1,'alphabet','nt','scoringmatrix',scmat};
maparams = {'gapopen',32,'extendgap',1,'scoringmatrix',scmat,'TERMINALGAPADJUST',true,...
 'weights','equal','delaycutoff',1e-10};

ths = [50 80 90 98];
fprintf('Aligning to amplicons:\n');
if isempty(P.range_to_process)
  nfirst=1; nlast=nr;
else
  nfirst=P.range_to_process(1);
  nlast=P.range_to_process(2);
end

fout = fopen([P.outdir '/output.txt'],'wt');

for r=nfirst:nlast;

  idx=find(X.amp==r); amp=R.amplicon{r}; if 1; sc_f = nan(length(idx),1); sc_r = sc_f; fprintf('\nr %d   %d seqs',r,length(idx));
  maxsc = sum(amp~='N'); for i=1:length(idx)
    [af{i} bf{i} cf{i}] = swalign(X.seq{idx(i)},amp,swparams{:}); sc_f(i) = 100*sum(bf{i}(2,:)=='|')/maxsc;
    [ar{i} br{i} cr{i}] = swalign(rc(X.seq{idx(i)}),amp,swparams{:}); sc_r(i) = 100*sum(br{i}(2,:)=='|')/maxsc;
  end, s=max(sc_f,sc_r);fprintf('   %d >50   %d >80   %d >90   %d >98\n',sum(s>50),sum(s>80),sum(s>90),sum(s>98)); end
  fprintf(fout,'\nr %d   %d seqs   %d >50   %d >80   %d >90   %d >98\n',r,length(idx),sum(s>50),sum(s>80),sum(s>90),sum(s>98));

  for thi=1:length(ths), th=ths(thi);
   try
    lim1 = th;
    seqs = ['k';'k';amp;'k';'k';X.seq(idx(sc_f>=lim1));'k';'k';rc(X.seq(idx(sc_r>=lim1)));'k';'k';amp;'k';'k'];
    scores = [inf(5,1);sc_f(sc_f>=lim1);inf(2,1);sc_r(sc_r>=lim1);inf(5,1)];
    regions =[nan(5,1);X.region(idx(sc_f>=lim1));nan(2,1);X.region(idx(sc_r>=lim1));nan(5,1)];
    tic;fprintf('multialigning %d seqs... ',length(seqs)-12); ma = multialign(seqs,maparams{:});toc
  
    lim2 = th;
    map('k-ACGTN')=[0 0 1:5]; ii=(scores>=lim2);q=map(ma(ii,:));jj=(regions(ii)-1==(r>192));q(jj,:)=q(jj,:)+6;
    % A=red C=blue G=green T=yellow N=white
    colormap(0.1*[0 0 0;6 2 2;2 2 6;2 6 2;6 6 2;10 10 10;1 1 2;7 3 3;3 3 7;3 7 3;7 7 3;10 10 10]);
    q(end,end)=11;imagesc(q);set(gca,'position',[0.03 0.05 0.94 0.9],'visible','off');

    print_to_file([P.outdir '/stringency' num2str(th) '_alignment' num2str(r) '.png']);
   catch me
    fprintf('Error: %s\n',me.message); 
   end
  end

  idx_f = idx(sc_f>=90); idx_r = idx(sc_r>=90);
  fprintf('%d    T/N  %d/%d      Tot  %d      Reg1  %d/%d     Reg2  %d/%d\n',r,...
    R.tumreads(r),R.normreads(r),length(idx),...
    sum(X.region(idx_f)==1),sum(X.region(idx_r)==1),...
    sum(X.region(idx_f)==2),sum(X.region(idx_r)==2));
  fprintf(fout,'%d    T/N  %d/%d      Tot  %d      Reg1  %d/%d     Reg2  %d/%d\n',r,...
    R.tumreads(r),R.normreads(r),length(idx),...
    sum(X.region(idx_f)==1),sum(X.region(idx_r)==1),...
    sum(X.region(idx_f)==2),sum(X.region(idx_r)==2));

  X.hit(idx_f) = 1; X.hit(idx_r) = -1;
end

fclose(fout);

% judge each rearrangement

rmin = min(X.region);
rmax = max(X.region);
Fhits = zeros(nr,rmax);
Rhits = zeros(nr,rmax);
for i=rmin:rmax
  Fhits(:,i) = histc(X.amp(X.hit==1 & X.region==i),1:nr);
  Rhits(:,i) = histc(X.amp(X.hit==-1 & X.region==i),1:nr);
end

R.tumFhits = [Fhits(1:192,2);Fhits(193:384,1)];
R.tumRhits = [Rhits(1:192,2);Rhits(193:384,1)];
R.normFhits = [Fhits(1:192,1);Fhits(193:384,2)];
R.normRhits = [Rhits(1:192,1);Rhits(193:384,2)];
R.tumhits = R.tumFhits + R.tumRhits;
R.normhits = R.normFhits + R.normRhits;

R.validation_result = repmat({'not_detected'},nr,1);
idx = find(R.tumhits>0);
R.validation_result(idx) = repmat({'somatic'},length(idx),1);
idx = find(R.normhits>0 & (R.tumhits./R.normhits)<12);
R.validation_result(idx) = repmat({'germline'},length(idx),1);




