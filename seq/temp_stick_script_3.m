function h = temp_stick_script_3(M,G,g,stickdir)
% temp_stick_script_3(M,G,g,stickdir)
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
% used in /cga/tcga-gsc/home/lawrence/mut/analysis/20130306_pancan12/run.m

gname = G.gene{g};

% load domain info from stickfigs feature file
txtname = [stickdir '/' gname '.stick_fig.txt'];
if ~exist(txtname,'file'), error('no features file for %s',gname); end % doesn't work without this file
f = load_struct(txtname);
f = make_numeric(f,{'dnapos','start','end','R','G','B'});

% remove mutations and domains with unknown aa_pos
f = reorder_struct_exclude(f,isnan(f.start)|isnan(f.end));

% if no "protein" line, supply one
idx=find(strcmp('protein',f.feature));
if isempty(idx)
  fprintf('No "protein" line: supplying one.\n');
  p=[]; p.feature = {'protein'};
  p.dnapos = 0; p.start = 1; p.end = max(f.end);
  p.label = {gname};
  f = concat_structs_keep_all_fields({f,p});
end

% get protein length
protlength = f.end(find(strcmp('protein',f.feature)));

% some edits to header
fh = reorder_struct_exclude(f,strcmp('mutation',f.feature)); fm = reorder_struct(f,strcmp('mutation',f.feature));
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

% extract original mutation information from M struct
mutidx = find(strcmp(M.mut.gene,gname));
m = reorder_struct(M.mut,mutidx);
m.mutidx = mutidx;
m.midx = (1:slength(m))';
% remove noncoding mutations
m = reorder_struct_exclude(m,grepmi('intron|igr|utr|flank',m.type));

% make unified mutation list
m.fmidx = multimap(m,fm,{'start','newbase'},{'dnapos','newbase'});
fm = rename_fields(fm,{'start','end','type','newbase'},{'aa_start','aa_end','type2','newbase2'});
m = merge_structs({m,reorder_struct(fm,m.fmidx)}); clear fm

% again, remove mutations and domains with unknown aa_pos
m = reorder_struct_exclude(m,isnan(m.aa_start)|isnan(m.aa_end));

% change shape of synonymous mutations to diamond
idx = find(grepmi('silent|synon',m.type2));
m.shape(idx) = repmat({'diamond'},length(idx),1);

% saturate missenses according to conservation (highly saturated = highly conserved)
is_mis = grepmi('missense',m.type2);
if any(is_mis)
  m = make_numeric(m,{'cons46'});
  sat=max(0,m.cons46-50)/50;
  m.R(m.R==0 & is_mis)=(1-sat(m.R==0 & is_mis));
  m.G(m.G==0 & is_mis)=(1-sat(m.G==0 & is_mis));
  m.B(m.B==0 & is_mis)=(1-sat(m.B==0 & is_mis));
end




% y-positions for display: segregate by tumor type
rand('twister',1234);
m.ypos = nan(slength(m),1);
if isfield(m,'ttype')
  % sort by tumor type
  [ttu ttui ttuj] = unique(m.ttype);
  divheight = 0.025;
  while(true)
    mutheight = (1-((length(ttu)-1)*divheight))/slength(m);
    if mutheight>0
      break
    else
      divheight = divheight/2;
    end
  end
  ystart = 1;
  for tti=1:length(ttu)
    mi = find(ttuj==tti);
    nmi = length(mi);
    ord = randperm(nmi);
    if 0&& nmi>2 %&& strcmp(ttu{tti},'Colorectal')    % DOESNT WORK YET
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % do a small optimization to reduce visual overlap of symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      q=[]; q.y=ord'; q.x=(m.aa_start(mi)/max(m.aa_start))*20*slength(m);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % randomly pairs trios of mutations,
      % swapping y1,y2 if that increases total separation
      fprintf('[Optimize %s:',ttu{tti});
      ncycles = min(10000,2*nmi^2);
      m1=1; m2=1;
      nswaps = 0;
      for cycle=1:ncycles
 %       if ~mod(cycle,1e4), fprintf(' %d/%d (%d)',ceil(cycle/1e5),ceil(ncycles/1e5),nswaps); nswaps=0; end
        m1=m1+1; if m1>nmi, m1=1; m2=m2+1; if m2>nmi, m2=1; end; end
        XY = [q.x q.y];
        thresh = 20;
        dist1 = sqrt(sum(bsxfun(@minus,XY,XY(m1,:)).^2,2));
        dist2 = sqrt(sum(bsxfun(@minus,XY,XY(m2,:)).^2,2));
        totdist_now = sum(dist1<thresh)+sum(dist2<thresh);
        tmp = XY(m2,2); XY(m2,2)=XY(m1,2); XY(m1,2)=tmp;
        dist1_swap = sqrt(sum(bsxfun(@minus,XY,XY(m1,:)).^2,2));
        dist2_swap = sqrt(sum(bsxfun(@minus,XY,XY(m2,:)).^2,2));
        totdist_swap = sum(dist1_swap<thresh)+sum(dist2_swap<thresh);
        if 00 & totdist_swap > totdist_now
          col = 0.8*ones(nmi,3);
          thresh = 12;
          idx = find(dist1<thresh); col(idx,:)=repmat([0.5 0 0],length(idx),1);
          idx = find(dist2<thresh); col(idx,:)=repmat([0 0 0.5],length(idx),1);
          col(m1,:) = [1 0 0]; % red
          col(m2,:) = [0 0 1]; % blue
          figure(1);clf,scatter(q.x,q.y+0.5,100,col);set(gca,'position',[0.05 0.4 0.9 0.5]);set(gcf,'position',[50 350 50+2.1*(max(q.x)/20) 350])
          fprintf('WILL SWAP> ');
          keyboard
        end

        if totdist_swap > totdist_now
          tmp = q.y(m1);
          q.y(m1) = q.y(m2);
          q.y(m2) = tmp;
          nswaps = nswaps + 1;
        end
      end
      fprintf(' done]\n');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ord = q.y';
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    yend = ystart - (nmi-1)*mutheight;
    m.ypos(mi(ord)) = (ystart:-mutheight:yend) + mutheight*0.5 - divheight*0.05;
    ystart = yend - mutheight - divheight;
  end
else
  % just sort randomly
  ord = randperm(slength(m));
  m.ypos(ord) = (1:slength(m))'/slength(m);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resave feature file
fm = reorder_struct(m,~isnan(m.fmidx));
fm = rename_fields(fm,{'aa_start','aa_end','type2','newbase2'},{'start','end','type','newbase'});
mut_add_fields = {'ttype','patient','chr','i_tumor_f','ref_allele'};
fm = keep_fields(fm,[fieldnames(fh);'ypos';as_column(mut_add_fields)]);
f = concat_structs_keep_all_fields({fh,fm});
txtname2 = [stickdir '/' gname '.stick_fig.3.txt']; save_struct(f,txtname2);


return


% other stuff to add later:

% DRAW SPLICESITES
if isfield(M,'gene') && isfield(M.gene,'intron_starts') && isfield(M.gene,'intron_ends')
  % find aa positions of splice sites
  mgidx = find(strcmp(M.gene.name,gname));
  if ~isempty(mgidx)
    ss = []; ss.dnapos = union(M.gene.intron_starts{mgidx},M.gene.intron_ends{mgidx});
    if ~isempty(ss.dnapos)
      ss.aapos = nan(slength(ss),1);
      % for each one, find the closest mutation and map according to that
      map = keep_fields(f,{'dnapos','start'}); map=rename_field(map,'start','aapos');
      map = reorder_struct_exclude(map,map.dnapos==0 | isnan(map.dnapos) | map.aapos==0 | isnan(map.aapos));
      map = sort_struct(map,'dnapos');
      for i=1:slength(ss)
        dist = ss.dnapos(i)-map.dnapos;
        [tmp j] = min(abs(dist));
        mapped_aapos = map.aapos(j)+dist(j);
        if (mapped_aapos>=1 && mapped_aapos<=protlength && dist(j)<=1000)
          ss.aapos(i) = mapped_aapos;
        end
      end
      % draw splice sites
      for i=1:slength(ss)
        if ~isnan(ss.aapos(i))
          x = fig.prot.left + (fig.prot.xscale.*(ss.aapos(i)-1));
          line(x*[1 1],[fig.prot.bottom fig.prot.top],'color',[1 0 0]);
        end
      end
    end
  end
end

% DRAW LINES TO DIVIDE TUMOR TYPES
if isfield(m,'ttype')
  % ADD:
  % -- lines to divide tumor types
  % -- tumor type names
  % -- copy number icons
  cn_available = (isfield(M,'cn') && isfield(M.cn,'pat') && isfield(M.cn.pat,'ttype'));
  ystart = 1;
  for tti=1:length(ttu)
    mi = find(ttuj==tti);
    yend = ystart - (length(mi)-1)*mutheight;
    if tti<length(ttu)
      ly = fig.mut.bottom + yend*fig.mut.height;
      line([fig.mut.left fig.mut.right],[ly ly],'color',[0 0 0]);
    end
    params = {'verticalalignment','middle'};
    if P.mutfig_style==1 || P.mutfig_style==2
      tx = fig.mut.left/2;
      ty = fig.mut.bottom + ((ystart+yend)/2+divheight/2)*fig.mut.height;
      params = [params,'horizontalalignment','center','fontsize',10];      
    elseif P.mutfig_style==3
      ty = fig.mut.bottom + ((ystart+yend)/2+0.85*divheight)*fig.mut.height-0.005;
      if cn_available
        tx = fig.prot.left * 0.75;
        params = [params,'horizontalalignment','right','fontsize',14];
      else
        tx = fig.prot.left * 0.95;
        params = [params,'horizontalalignment','right','fontsize',18];
      end
    end
    if isfield(M.ttype,'color')
      ttidx = find(strcmp(ttu{tti},M.ttype.name));
      if ~isempty(ttidx)
        params = [params,'color',M.ttype.color(ttidx,:)];
      end
    end
    text(tx,ty,ttu{tti},params{:});
    if isfield(M,'cn') && isfield(M.cn,'pat') && isfield(M.cn.pat,'ttype')
      % DRAW COPY-NUMBER BAR
      gidx = find(strcmp(gname,M.cn.gene.name));
      pidx = find(strcmp(ttu{tti},M.cn.pat.ttype));
      if length(gidx)==1 && ~isempty(pidx)
        cn = sort(M.cn.dat(gidx,pidx));
        % expect to receive copy number data on scale from about -1 to +1
        cn(cn<-1) = -1;
        cn(cn>+1) = +1;
        cn_idx = floor((cn+1)*32)+1;
        cn_idx(cn_idx<1) = 1;
        cn_idx(cn_idx>64) = 64;
        bwr_colors = bwr();
        % draw CN bar
        fig.cn.left = fig.prot.left * 0.77;
        fig.cn.right = fig.prot.left * 0.98;
        fig.cn.width = fig.cn.right-fig.cn.left;
        fig.cn.height = 0.02;
        top = ty + fig.cn.height/2;
        bottom = ty - fig.cn.height/2;
        rectangle('position',[fig.cn.left bottom fig.cn.width fig.cn.height],'edgecolor',[0 0 0]);
        for i=1:length(cn)
          rectangle('position',[fig.cn.left+fig.cn.width*((i-1)/length(cn)) bottom fig.cn.width/length(cn) fig.cn.height],...
                  'linestyle','none','facecolor',bwr_colors(cn_idx(i),:));
        end  
      end
    end
    ystart = yend - mutheight - divheight;
  end
end




