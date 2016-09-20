function [m use_midx control_midx] = temp_stick_script(M,G,g)
% temp_stick_script(M,G,g)
%
% wrapper for mutfig.m
% --adds new visualization elements
%
% M = master mutation data object
% G = genelist
% g = which gene in G to work with
%
% e.g.
% load('/cga/tcga-gsc/home/lawrence/mut/analysis/20121031_pancan/data.v2.2.mat','M');
% load('/cga/tcga-gsc/home/lawrence/mut/analysis/20121031_pancan/genes.1.mat','G','ttypes');
%
% used in /cga/tcga-gsc/home/lawrence/mut/analysis/20121031_pancan/run.m

use_midx = [];
control_midx = [];

% pull out mutations of our gene
%look(G,g);
gname = G.gene{g};
mutidx = find(strcmp(M.mut.gene,gname));
m = reorder_struct(M.mut,mutidx);
m.mutidx = mutidx;
m.midx = (1:slength(m))';

% load domain info from stickfigs feature file
stickdir = '/cga/tcga-gsc/home/lawrence/mut/analysis/20121031_pancan/stickfigs/run1';
txtname = [stickdir '/' gname '.stick_fig.txt'];
if ~exist(txtname,'file'), return; end % doesn't work without this file
f = load_struct(txtname,'%s%f%s%s%f%f%s%f%f%f%s%s');
fh = reorder_struct_exclude(f,strcmp('mutation',f.feature)); fm = reorder_struct(f,strcmp('mutation',f.feature));

% some edits to header
idx=find(strcmp('protein',fh.feature)); fh.R(idx)=0.8; fh.G(idx)=0.8; fh.B(idx)=0.8;
aamax = max(f.end); fh.end(idx) = aamax;

% for domains with same label, color them all the same color
d=find(strcmp('domain',fh.feature));
if ~isempty(d)
  [u ui uj] = unique(fh.label(d));
  for j=1:length(u)
    idx = d(uj==j);
    fh.R(idx) = fh.R(idx(1));
    fh.G(idx) = fh.G(idx(1));
    fh.B(idx) = fh.B(idx(1));
  end
end

% make unified mutation list
m.fmidx = multimap(m,fm,{'start','newbase'},{'dnapos','newbase'});
fm = rename_fields(fm,{'start','end','type','newbase'},{'aa_start','aa_end','type2','newbase2'});
m = merge_structs({m,reorder_struct(fm,m.fmidx)}); clear fm

% make fix for weird cases: a few frameshifts that are annotated like missenses?
idx = grep('fs',m.Protein_Change,1); 
m.is_null(idx)=1; m.is_missense(idx)=0;
m.is_silent(idx)=0; m.is_inframe(idx)=0; m.is_indel(idx)=1;

% saturate according to conservation (highly saturated = highly conserved)
sat=max(0,m.cons46-50)/50;
m.R(m.R==0)=(1-sat(m.R==0));
m.G(m.G==0)=(1-sat(m.G==0));
m.B(m.B==0)=(1-sat(m.B==0));

% prioritize mutations
m.prior = m.mutdens .* m.cons46;
m.prior(m.is_null>0) = 0;
m.prior(m.is_silent>0 | m.is_inframe>0) = -1;

% y-positions for display
[tmp ord] = sort(m.prior);
m.ypos = nan(slength(m),1);
m.ypos(ord) = (1:slength(m))'/slength(m);

% resave feature file and draw figure
fm = reorder_struct(m,~isnan(m.fmidx));
fm = rename_fields(fm,{'aa_start','aa_end','type2','newbase2'},{'start','end','type','newbase'});
fm = keep_fields(fm,[fieldnames(fh);'ypos']);
f = concat_structs_keep_all_fields({fh,fm});
txtname2 = [stickdir '/' gname '.stick_fig.2.txt']; save_struct(f,txtname2);
P=[]; P.mutfig_style = 2; clf, fig = draw_mutfig(txtname2,P);

%%% draw POST-MUTFIG ADDITIONS
%%% AND pull out the mutations we want to test

m2 = reorder_struct(m,~isnan(m.aa_start));
m2 = sort_struct(m2,'aa_start');
m2 = reorder_struct(m2,m2.is_missense>0); % missenses only
m2 = reorder_struct(m2,m2.start==m2.end); % SNPs only

if isfield(G,'cmut')
  nmut = G.cmut(g);
else
  nmut = 8;
end

mutno=0;
use_midx = [];
finished = false;
while (~finished)

  % make track that is the smoothed distribution of sum-of-conservation
  h = zeros(aamax,1);
  for i=1:slength(m2)
    hi = m2.aa_start(i);
    h(hi) = h(hi) + m2.cons46(i);
  end
  smoothwin = ceil(aamax/20);
  h = smooth2(h,smoothwin);
  hmax = max(h);
  hold on
  %  plot(fig.mut.left + fig.mut.xscale * (1:aamax), fig.mut.bottom + fig.mut.height * (h/hmax),'-b','linewidth',2);
  hold off

  % find highest-priority mutation
  mutno = mutno + 1;
  [tmp idx] = max(m2.prior);
  midx = m2.midx(idx);
  use_midx = [use_midx; midx];

  % reduce priority in a window around the chosen mutation
  reduce_radius = 10;
  reduce_factor = 2;
  dist = abs(m2.aa_start - m2.aa_start(idx));
  rf = reduce_factor * min(1,(reduce_radius./dist));
  fac = min(1,1./rf);
  m2.prior = m2.prior .* fac;

  % completely remove from further consideration all mutations that have the same Protein_Change label
  m2 = reorder_struct_exclude(m2,strcmp(m2.label,m2.label{idx}));

  if (mutno==nmut), finished = true; end
end

% choose the control mutation

if nmut>=3
  m3 = reorder_struct(m,~isnan(m.aa_start));
  m3 = sort_struct(m3,'aa_start');
  m3 = reorder_struct(m3,m3.is_missense>0); % missenses only
  m3 = reorder_struct(m3,m3.start==m3.end); % SNPs only

  h = histc(m3.aa_start,1:aamax);   % histogram of *exact* positions

  m3 = reorder_struct_exclude(m3,ismember(m3.aa_start,m.aa_start(use_midx)));  % remove positions that were selected for the experimental mutations

  % priority here is (100-conservation), divided by number of times (squared) that exact mutation occurs
  m3.ntimes = nansub(h,m3.aa_start);
  m3.prior = (100-m3.cons46) ./ (m3.ntimes.^2);
  % hard threshold: make sure conservation is <60
  consthresh=60;
  m3 = reorder_struct(m3,m3.cons46<consthresh);
  if slength(m3)>0
    % add a bonus for proximity to muts 1 and 2
    dist1 = abs(m3.aa_start-(m.aa_start(use_midx(1))));
    dist2 = abs(m3.aa_start-(m.aa_start(use_midx(2))));
    m3.dist = min(dist1,dist2);
    m3.pd = m3.prior + 10*(m3.dist<(aamax/5)) + 10*(m3.dist<(aamax/3));
    m3 = sort_struct(m3,'pd',-1);
    control_midx = m3.midx(1);
    % remove one of the experimental mutations to make room for the control mutation
    use_midx(end)=[];
  else
    control_midx = [];
  end
end

%%%%%%%%% DISPLAY THE LOCATIONS OF THE CHOSEN MUTATIONS

if ~isempty(use_midx)
  % display the experimental mutations
  hold on
  for mutno=1:length(use_midx), midx=use_midx(mutno);
    x = fig.mut.left + fig.mut.xscale * (m.aa_start(midx)-1);
    y = fig.mut.top + 0.05 + 0.08 * ((mutno-1)/nmut);
    plot(x,y, '.g','markersize',30);
    text(x,y,num2str(mutno),'horizontalalignment','center');
  end
  hold off
  pr(m,{'gene','label','cons46'},use_midx);
end

if ~isempty(control_midx)
  disp('-------------------');
  % display the control mutation
  midx = control_midx;
  hold on
  x = fig.mut.left + fig.mut.xscale * (m.aa_start(midx)-1);
  y = fig.mut.top + 0.05 + 0.08 + 0.01;
  plot(x,y, '.','color',[0.8 0.8 0.8],'markersize',30);  % grey
  text(x,y,'c','horizontalalignment','center');
  hold off
  pr(m,{'gene','label','cons46'},control_midx);
end






