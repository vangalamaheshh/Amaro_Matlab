function [res M] = radial_spectrum_plot(M,P)
% radial_spectrum_plot(M,P)
% 
% based on mutation_spectra_pca_plot.m
%
% Mike Lawrence 2011-2013

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'build','hg19');
P = impose_default_value(P,'categdir',['/cga/tcga-gsc/home/lawrence/db/' P.build '/context65']);
P = impose_default_value(P,'impute_full_coverage',false);
P = impose_default_value(P,'pca_plot_output_filestem',[]);
P = impose_default_value(P,'collapse_strands_first',true);
P = impose_default_value(P,'method','NMF');
P = impose_default_value(P,'nfactors',6);
P = impose_default_value(P,'randseed',1240);
P = impose_default_value(P,'data','relrate');
P = impose_default_value(P,'manually_collapse_factors',[]);
P = impose_default_value(P,'radial_plot_by_clustering',false);
P = impose_default_value(P,'scale_NMF_factors_and_weightings',true);
P = impose_default_value(P,'find_names_of_categories',true);
P = impose_default_value(P,'max_k_for_category_name_finding',4);
P = impose_default_value(P,'optimize_for','screen'); % screen or print

% minimal required input data
% see pr/analysis/20111112_pca/run.m for example of running from a FH MutSig output

if ischar(M) || (isstruct(M) && ~isfield(M,'mut') && isfield(M,'gene'))
  % this is a maf
  M = maftoM(M);
end

if isfield(M,'patient') && ~isfield(M,'pat'), M=rename_field(M,'patient','pat'); end
if isfield(M.mut,'patient_idx') && ~isfield(M.mut,'pat_idx'), M.mut=rename_field(M,'patient_idx','pat_idx'); end
demand_fields(M,{'pat','mut'});if ~isfield(M.pat,'name') && isfield(M.pat,'indiv'), M.pat.name = M.pat.indiv; end
demand_fields(M.pat,{'name','ttype'});
demand_fields(M.mut,{'newbase'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(M.mut,'context65')
  M.mut.context65 = get_context(M.mut.chr,M.mut.start,P.categdir,P);
else
  M.mut = make_numeric(M.mut,'context65');
end

if ~isfield(M,'ttype'), M.ttype.name = unique(M.pat.ttype); end
if ~isfield(M.pat,'ttype_idx')
  M.pat.ttype_idx = listmap(M.pat.ttype,M.ttype.name);
else
  M.pat = make_numeric(M.pat,'ttype_idx');
end
if ~isfield(M.mut,'pat_idx')
  M.mut.pat_idx = listmap(M.mut.patient,M.pat.name);
else
  M.mut = make_numeric(M.mut,'pat_idx');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = generate_categ_context65_names();

if ~isfield(M,'Nn') 
  fprintf('Preparing data...\n');

  M.np = slength(M.pat);
  
  if isfield(M.mut,'classification'), M.mut = reorder_struct(M.mut,strcmp('SNP',M.mut.classification)); end
  M.mut.newbase_idx = listmap(M.mut.newbase,{'A','C','G','T'});
  M.mut = reorder_struct(M.mut,~isnan(M.mut.newbase_idx));
  if slength(M.mut)==0, error('No mutations left after filtering on M.mut.classification and M.mut.newbase'); end

  if isfield(M,'gene')
    M.ng = slength(M.gene);
    if ~isfield(M.mut,'gene_idx')
      if isfield(M.gene,'name') && isfield(M.mut,'gene')
        M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);
      end
    end
  end
  if isfield(M.mut,'gene_idx') && isfield(M.gene,'strand')
    M.mut = make_numeric(M.mut,'gene_idx');
    M.mut = reorder_struct(M.mut,~isnan(M.mut.gene_idx));
    cidx = listmap(M.gene.name,M.cov.gene.name);
    plusstrand = strcmp(M.gene.strand,'+');
    minusstrand = strcmp(M.gene.strand,'-');
    M.mut.is_plus=false(slength(M.mut),1);
    M.mut.is_plus(plusstrand(M.mut.gene_idx))=true;
    M.mut.is_minus=false(slength(M.mut),1);
    M.mut.is_minus(minusstrand(M.mut.gene_idx))=true;
    strand_info_available = true;
  else
    fprintf('No gene/strand information available: will not be able to orient to sense strand\n');
    strand_info_available = false;
  end

  categs_txt = [P.categdir '/categs.txt'];

  compbase('ACGT') = 'TGCA';
  idx1=nan(64,1); idx2=idx1;
  for i=1:64
    oldname = X.name{i};
    newname = [compbase(oldname(1)) ' in ' compbase(oldname(end)) '_' compbase(oldname(end-2))];
    idx1(i) = find(strcmp(oldname,X.name));
    idx2(i) = find(strcmp(newname,X.name));
  end
  
  M.Nn = nan(64,5,M.np);

  if ~P.impute_full_coverage
    if ~isfield(M,'cov') || ~isfield(M.cov,'gene_samp_orig_cov')
      fprintf('Actual coverage not available: will impute full coverage\n');
      P.impute_full_coverage = true;
    end
  end

  if ~P.impute_full_coverage
    for si=1:M.np, if ~mod(i,100), fprintf('%d/%d ',si,M.np); end
      fprintf('Please double-check correctness of this code...\n');
      keyboard
      %%%%%%% NEED TO CAREFULLY CHECK THIS
      cov = collapse_categories_to_65(M.cov.gene_samp_orig_cov,categs_txt);
%      N = squeeze(M.cov(:,si,1:64))';
      N = squeeze(cov(:,si,1:64))';
      Nplus = sum(N(:,plusstrand),2);
      Nminus = sum(N(:,minusstrand),2);
      idx = find(M.mut.pat_idx==si);
      is_pat = (M.mut.pat_idx==si);
      plusidx = (is_pat & M.mut.is_plus);
      minusidx = (is_pat & M.mut.is_minus);
      nplus = hist2d_fast(M.mut.context65(plusidx),M.mut.newbase_idx(plusidx),1,64,1,4);
      nminus = hist2d_fast(M.mut.context65(minusidx),M.mut.newbase_idx(minusidx),1,64,1,4);
      M.Nn(:,1,si) = Nplus(idx1) + Nminus(idx2);
      M.Nn(:,2,si) = nplus(idx1,1) + nminus(idx2,4);
      M.Nn(:,3,si) = nplus(idx1,2) + nminus(idx2,3);
      M.Nn(:,4,si) = nplus(idx1,3) + nminus(idx2,2);
      M.Nn(:,5,si) = nplus(idx1,4) + nminus(idx2,1);
    end, if M.np>100, fprintf('\n'); end
    
    %    n = zeros(64,4,M.ng); idx = find(M.mut.pat_idx==si);
    %    for k=1:length(idx), j=idx(k);
    %      a=M.mut.context(j); b = M.mut.newbase_idx(j); c = M.mut.gene_idx(j);
    %      n(a,b,c) = n(a,b,c) + 1;
    %    end
    %    Nngc = zeros(64,5,M.ng); Nngc(:,1,:) = N; Nngc(:,2:5,:) = n;
    %    M.Nn(:,:,si) = collapse_genes_to_coding_strand(Nngc,M.gene);
    %  end,fprintf('\n');

  else   % impute full coverage

    if isfield(M,'cov') && isfield(M.cov,'gene_orig_terr')
%       (low-memory-requirement option:)
%       terr = collapse_categories_to_65_tr(M.cov.gene_orig_terr,categs_txt);
      N = collapse_categories_to_65(M.cov.gene_orig_terr',categs_txt);
      N = N(1:64,:);
    else
      fprintf('No coverage information available: will use average exome coverage.\n');
      N = average_exome_coverage();
      strand_info_available = false;
    end

    if strand_info_available
      % reorient all mutations and coverage with respect to the sense strand
      Nplus = sum(N(:,plusstrand),2);
      Nminus = sum(N(:,minusstrand),2);
      for si=1:M.np, if ~mod(si,100), fprintf('%d/%d ',si,M.np); end
        is_pat = (M.mut.pat_idx==si);
        plusidx = (is_pat & M.mut.is_plus);
        minusidx = (is_pat & M.mut.is_minus);
        nplus = hist2d_fast(M.mut.context65(plusidx),M.mut.newbase_idx(plusidx),1,64,1,4);
        nminus = hist2d_fast(M.mut.context65(minusidx),M.mut.newbase_idx(minusidx),1,64,1,4);
        M.Nn(:,1,si) = Nplus(idx1) + Nminus(idx2);
        M.Nn(:,2,si) = nplus(idx1,1) + nminus(idx2,4);
        M.Nn(:,3,si) = nplus(idx1,2) + nminus(idx2,3);
        M.Nn(:,4,si) = nplus(idx1,3) + nminus(idx2,2);
        M.Nn(:,5,si) = nplus(idx1,4) + nminus(idx2,1);
      end, if M.np>100, fprintf('\n'); end

    else  % no strand info available: can't reorient to sense strand
      Nall = sum(N,2);
      for si=1:M.np, if ~mod(si,100), fprintf('%d/%d ',si,M.np); end
        is_pat = (M.mut.pat_idx==si);
        n = hist2d_fast(M.mut.context65(is_pat),M.mut.newbase_idx(is_pat),1,64,1,4);
        M.Nn(:,1,si) = Nall(idx1);
        M.Nn(:,2:5,si) = n(idx1,:);
      end, if M.np>100, fprintf('\n'); end
    end

  end

  if P.collapse_strands_first
    M.Nn = collapse_Nn_64_by_strand(M.Nn);
  end

end

% remove samples with extremely few (<10) mutations
M.pat.Ntot = squeeze(sum(sum(M.Nn(:,1,:),1),2));
M.pat.ntot = squeeze(sum(sum(M.Nn(:,2:5,:),1),2));
idx = find(M.pat.ntot>=10);
M.pat.orig_idx = (1:slength(M.pat))';
M.pat = reorder_struct(M.pat,idx);
M.Nn = M.Nn(:,:,idx);
M.np = slength(M.pat);
M.mut.pat_idx = mapacross(M.mut.pat_idx,M.pat.orig_idx,(1:slength(M.pat))');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data

res=[];
res.pat = M.pat;
res.ttype = M.ttype;

[res.Nn res.cat] = toggle_Nn_5col_2col(M.Nn);
N = squeeze(res.Nn(:,1,:));
n = squeeze(res.Nn(:,2,:));

mutfrac = bsxfun(@rdivide,n,sum(n,1));
rate = n./N;
totrate = sum(n,1)./sum(N,1);
relrate = bsxfun(@rdivide,rate,totrate);
lograte = min(6,-log10(rate));
logrelrate = -log10(relrate); logrelrate(isinf(logrelrate)) = mean(logrelrate(~isinf(logrelrate)));

% softmax
softmax = -log10(relrate);
softmax(isinf(softmax)) = nan;
samplemean = nanmean(softmax,1);
softmax(isnan(softmax)) = 0;
softmax = bsxfun(@minus,softmax,samplemean);

switch P.data
    case 'count', data=n;
    case 'frac', data=mutfrac;
    case 'rate', data=rate;
    case 'relrate', data=relrate;
    case 'lograte', data=lograte;
    case 'logrelrate', data=logrelrate;
    case 'softmax', data=softmax;
    otherwise, error('invalid P.data = %s',P.data);
end
res.pat = M.pat;
res.data = data';

switch P.method
    case 'PCA' 
        rand('twister',P.randseed);
        c = princomp(data');
        %c(:,1) = -c(:,1);  % (to flip "USA" from softmax left-to-right)
        c(:,2) = -c(:,2);  % (to flip "USA" from softmax upsidedown)
        res.c = c;
        proj = c'*data;
    case 'FA'
        [fa.lambda, fa.psi, fa.t, fa.stats, fa.f] = factoran(data', P.nfactors);
        rand('twister',P.randseed);
        c = princomp(fa.f);
        res.fa = fa;
        res.c = c;
        proj = c'*fa.f';
    case 'NMF'
        rand('twister',P.randseed);
        [w,h]=nmf(data',P.nfactors,1);
        if isfield(P,'manually_collapse_factors') && ~isempty(P.manually_collapse_factors)
          w0=w; h0=h; nf0=P.nfactors; nf=length(P.manually_collapse_factors);
          fprintf('Manual factor collaspe from %d to %d factors.\n',nf0,nf);
          % collapse h
          h = nan(nf,size(h0,2));
          for i=1:nf
            j = P.manually_collapse_factors{i}; % components to add
            h1 = h0(j(1),:);
            max1 = max(h1(:));
            h(i,:) = h1;
            for k=2:length(j)
              hk = h0(j(k),:);
              maxk = max(hk(:));
              % scale each addend after the first one to the max of the first one
              hk_scaled = hk * (max1/maxk);
              h(i,:) = h(i,:) + hk;
%              h(i,:) = h(i,:) + hk_scaled;
            end
          end
          % recompute w
          w = data'/h;
          w(w<0) = 0;
          % reassign nf
          P.nfactors = nf;
        end
        if P.scale_NMF_factors_and_weightings
          hscale=1./sum(h,2);
          h=h.*repmat(hscale,1,size(h,2));
          w=w./repmat(hscale',size(w,1),1);
          wn=w./repmat(sum(w,2),1,size(w,2));
          wn(wn<0.01)=0; 
          res.wn=wn;
        end
        res.w=w; res.h=h;
        c = princomp(w);
        proj = c'*w';
    otherwise, error('invalid P.method = %s',P.method);
end

if isfield(res,'wn'), w=res.wn; else w=res.w; end

% Analyze factors further

ttpat = false(slength(M.pat),slength(M.ttype));
for i=1:slength(M.pat)
  tti=M.pat.ttype_idx(i);
  if tti>=1 && tti<=slength(M.ttype)
    ttpat(i,M.pat.ttype_idx(i))=true;
  end
end
if ~isfield(M.ttype,'nmut'), M.ttype.nmut = histc(nansub(M.pat.ttype_idx,M.mut.pat_idx),1:slength(M.ttype)); end
if ~isfield(M.pat,'nmut'), M.pat.nmut = histc(M.mut.pat_idx,1:slength(M.pat)); end

res.patfrac = nan(slength(M.pat),P.nfactors);
res.ttypefrac = nan(slength(M.ttype),P.nfactors);
res.ttypefrac_text = cell(P.nfactors,1);
res.breakdown = cell(P.nfactors,1);

for fi=1:P.nfactors
  fprintf('\n\nAnalyzing factor %d/%d\n\n',fi,P.nfactors);

  % find which mutations are explained by this factor
  terr = squeeze(res.Nn(:,1,:));
  totterr = sum(terr,1)';
  totmut = squeeze(sum(res.Nn(:,2,:),1));
  totrate = totmut ./ totterr;
  relrate_fi = res.w(:,fi)*res.h(fi,:);
  rate_fi = bsxfun(@times,relrate_fi,totrate);
  mut_fi = rate_fi.*terr';
  totmut_fi = sum(mut_fi,2);
  patfrac_fi = max(0,totmut_fi./M.pat.nmut);
  totmut_ttype_fi = ttpat'*totmut_fi;
  ttypefrac_fi = max(0,totmut_ttype_fi./M.ttype.nmut);
  res.patfrac(:,fi) = patfrac_fi;
  res.ttypefrac(:,fi) = ttypefrac_fi;
  
  % create sorted text version of tumortype fraction-explained
  [tmp ord] = sort(ttypefrac_fi,'descend'); res.ttypefrac_text{fi} = cell(length(ord),1);
  for j=1:length(ord),i=ord(j);
    res.ttypefrac_text{fi}{j}=sprintf('%15s  (%.0f%% of all mutations)',M.ttype.name{i},100*ttypefrac_fi(i));
  end

  if P.find_names_of_categories
    % analyze categories
    Nn_fi = round([sum(terr,2) sum(mut_fi,1)']);
    Nn_fi = toggle_Nn_5col_2col(Nn_fi);
    if ~P.collapse_strands_first
      Nn_fi = collapse_Nn_64_by_strand(Nn_fi);
    end
    PP=[]; PP.max_k = P.max_k_for_category_name_finding;
    res.breakdown{fi} = find_mut_categs(Nn_fi,PP);
    
    % assign names to factors (by looking at where the maximum discrimination is achieved)
    rrr = nan(PP.max_k,1);
    for k=3:PP.max_k, rrr(k) = res.breakdown{fi}{k}.relrate(1)/res.breakdown{fi}{k}.relrate(2); end
    [tmp kidx] = max(rrr);
    res.names{fi,1} = res.breakdown{fi,1}{kidx}.name{1};
    fprintf('Name chosen:  %s\n',res.names{fi});
  else
    res.names{fi,1} = ['Factor ' num2str(fi)];
  end
  if isfield(P,'force_names') && ~isempty(P.force_names)
    idx = find(P.force_names{1}==fi);
    if ~isempty(idx)
      idx=idx(1);
      res.names{fi} = P.force_names{2}{idx};
      fprintf('Manually specified name:  %s\n',res.names{fi});
    end
  end
  res.names2{fi,1} = [res.names{fi} ' (' num2str(fi) ')'];
end % (next factor)

if isfield(res.ttype,'long'),ttname=res.ttype.long; else ttname=res.ttype.name; end

% RADIAL PLOT: BY SORTING FACTORS   (= Figure 2 of the heterogeneity paper)

nf = P.nfactors;

if strcmpi(P.optimize_for,'print')
  P = impose_default_value(P,'pointsize',5);
else % screen
  P = impose_default_value(P,'pointsize',20);
end
P = impose_default_value(P,'radial_plot_gap',0);
P = impose_default_value(P,'radial_plot_sep',0.03);
P = impose_default_value(P,'radial_plot_jitter',0.1);
P = impose_default_value(P,'radial_plot_improved_allocation',true);
P = impose_default_value(P,'radial_plot_labelangularpos',0.50);
P = impose_default_value(P,'nudge_points_for_nature_figure',false);

% rotate the plot so the mutations/Mb axis will be at 3 o'clock
% (this is very hacky...need to find an automatic way of computing this from nf)
if nf==4, P = impose_default_value(P,'radial_plot_rotate',0.1);
elseif nf==5, P = impose_default_value(P,'radial_plot_rotate',-0.30);
elseif nf==6, P = impose_default_value(P,'radial_plot_rotate',-0.95);
elseif nf==7, P = impose_default_value(P,'radial_plot_rotate',-1.25);
elseif nf==8, P = impose_default_value(P,'radial_plot_rotate',-1.48);
elseif nf==9, P = impose_default_value(P,'radial_plot_rotate',-1.75);
elseif nf==10, P = impose_default_value(P,'radial_plot_rotate',-1.99);
elseif nf==12, P = impose_default_value(P,'radial_plot_rotate',-2.20);
elseif nf==15, P = impose_default_value(P,'radial_plot_rotate',-2.41);
else, P = impose_default_value(P,'radial_plot_rotate',0);
end

gap=P.radial_plot_gap; sep=P.radial_plot_sep; jitter=P.radial_plot_jitter; rotate=P.radial_plot_rotate;
rotate = rotate + ((sep-0.20)/2);  % (so that rotate is independent of sep)

frac = res.patfrac;
[tmp factord] = sort(frac,2,'descend'); method=1;

if P.radial_plot_by_clustering
  method=4;
elseif P.radial_plot_improved_allocation
  method=3;
else
  method=1;
end

if method==3 && nf>10
  fprintf('nf=%d is too large for current implementation of improved angular allocation: will use older method\n',nf);
  method=1;
end

if method==1
  %% (was used in the initial submission of the heterogeneity paper)
  theta=rotate; for i=1:nf
    if i==1, angle=(1-gap); elseif i==2, angle=(1-sep); else angle=1; end
    theta = theta + ((factord(:,i)-1))*angle*2*pi/(nf^i);
  end
  % add jitter
  theta = theta + jitter*rand(size(theta)) - jitter/2;
elseif method==2
  %% (was not used)
  [u ui uj] = unique(factord,'rows'); theta = ((uj-1)/length(u))*2*pi;
  % add jitter
  theta = theta + jitter*rand(size(theta)) - jitter/2;
elseif method==3
  %% improved allocation that doesn't waste space due to leaving a spot for first_factor==second_factor, etc.
  % convert factord to tag
  mult = (nf.^((0:(nf-1))'));
  factord_tag = factord * mult;
  % make list of all unique orderings and their tags
  O = (nf+1)-perms(1:nf);
  no = size(O,1);
  O_tag = O * mult;
  % assign angles to each ordering
  O_theta = nan(no,1);
  tot_radians_avail = (2*pi)*(1-(gap+sep));
  radians_per_ordering = tot_radians_avail / no;
  angle = rotate+gap+(sep/4);
  for i=1:no
    O_theta(i) = angle;
    if mod(i,no/nf)==0, angle = angle + sep; end
    angle = angle + radians_per_ordering;
  end
  % map factord to angle
  theta = mapacross(factord_tag,O_tag,O_theta);
  % add jitter
  theta = theta + jitter*rand(size(theta)) - jitter/2;
  if P.nudge_points_for_nature_figure
    % nudge two points in 8 o'clock sector to save horizontal space
    idx = find(rho>=3.4 & theta>3.4 & theta<3.6); rho(idx) = rho(idx)-0.3;
  end
elseif method==4
  % "RADIAL PLOT #2", by clustering the factor weights
  frac = res.patfrac; frac = frac./repmat(sum(frac,2),1,size(frac,2));
  D = make_D(frac);
  addpath /cga/tcga-gsc/home/lawrence/cga/trunk/matlab/POS
  [res.Dord,res.dend] = one_way_clustering(D,'row',struct('cluster','average','dist','correlation'));
  [tmp res.pat.dendord]=sort(res.dend.idx');
  clear pi
  theta=((res.pat.dendord-1)/length(res.pat.dendord))*(2*pi);
end

rho = -0.6+max(0,min(4,log10(res.pat.ntot)));
[x y] = polar_to_cartesian(rho,theta);

% assign colors to tumor types
M.ttype.color = get_ttype_colors(M.ttype.name);
res.ttype = M.ttype; % (to save colors)

figure(2),clf,hold on
% first plot points by tumor type, so legend will have correct colors
for i=1:slength(res.ttype),idx=find(M.pat.ttype_idx==i);
  scatter(x(idx),y(idx),P.pointsize,res.ttype.color(i,:),'filled');
end
% then re-plot points in random order to show a random sampling of tumor types on top
ord = randperm(M.np);
for i=1:length(ord),idx=ord(i);
  scatter(x(idx),y(idx),P.pointsize,res.ttype.color(res.pat.ttype_idx(idx),:),'filled');
end
hold off,set(gca,'visible','off','position',[0.05 0.12 0.76 0.76]);
w = max(abs([min(x) max(x) min(y) max(y)]));
P = impose_default_value(P,'radial_plot_labelradialpos',w);
xlim([-w w]);ylim([-w w]);
legend(ttname,'location','eastoutside','fontsize',8), legend('boxoff');
%set(gcf,'color',[1 1 1],'position',[292 25 1130 793]);
set(gcf,'color',[1 1 1],'position',[456   190   765   502]);

if method~=4
  % display category labels
  labelangularpos=P.radial_plot_labelangularpos; labelradialpos=P.radial_plot_labelradialpos;  % category labels
  for i=1:nf
    theta1=rotate+(2*pi*((i-1)/nf))-(sep/2);
    [x0 y0]=polar_to_cartesian(0.3,theta1); [x1 y1]=polar_to_cartesian(3,theta1); 
    line([x0 x1],[y0 y1],'color',0.7*ones(1,3),'linewidth',1);
    [x2 y2]=polar_to_cartesian(labelradialpos,theta1+labelangularpos*(2*pi*(1-gap)/nf));
    name = res.names{i};
    text(x2,y2,name,'horizontalalignment','center','fontsize',12);
  end
  % and mutation rate axis
  theta0=rotate+(2*pi*((round(nf/2)-2)/nf))-(sep/2);
  for log10rate=-6:1:-4
    log10n = log10((10.^log10rate)*median(res.pat.Ntot));
    rho0 = -0.6+max(0,min(4,log10n)); [x0 y0]=polar_to_cartesian(rho0,theta0);
    rectangle('position',[x0-0.01 y0-0.03 0.02 0.06],'facecolor',0.7*ones(1,3),'linestyle','none');
    text(x0,y0+0.15,num2str((10.^log10rate)*1e6),'horizontalalignment','center','fontsize',8)
    if x0>0.8*w, break; end % out of display range
  end
  text(x0-0.29,y0-0.15,'mutations/Mb','horizontalalignment','center','fontsize',8);
end

%%%%%%%%% NOTE:  MAY NEED TO BE REPEATED TO GET PROPER ASPECT RATIO  (PROBLEM WITH XWin-32 v10)

% DISPLAY SUMMARY

fprintf('\nFactors\n');
pr(res.names2);

fprintf('\nTumor type breakdowns by factor\n');
tmp=[]; tmp.ttype = res.ttype.name;
for f=1:size(res.ttypefrac,2), tmp.(['f' num2str(f)])=round(100*res.ttypefrac(:,f)); end
pr(tmp)




function N = average_exome_coverage
  N=[...
      516645      336113      545347      365625      462172      488457      797766      489731
      549735      357085      540864      397040      230606      269162      226201      252346
      469526      429764      160802      389408      664422      477669      235217      616696
      509693      512648      186053      541099      537050      568802      169936      589050
      577447      538629      601391      389620      168530      186767      234416      162102
      556167      509431      469219      426880      536918      514744      661776      473758
      252140      397574      486994      367907      226652      552864      803001      556516
      267878      359998      486805      338852      231118      560177      467349      524231
    ];
  N=N(:);

