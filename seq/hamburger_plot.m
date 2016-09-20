function M = hamburger_plot(M,P)
% Gaddy Getz and Mike Lawrence

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'hamburger_plot_output_filename',[]);
P = impose_default_value(P,'hamburger_sample_sort_method','rate');
P = impose_default_value(P,'include_edits_for_print_version',true);

if ischar(M) || (~isfield(M,'mut') && isfield(M,'gene'))
  % this is a maf filename or loaded maf struct
  M = maftoM(M);
end

if isfield(M,'patient') && ~isfield(M,'pat'), M = rename_field(M,'patient','pat'); end
if isfield(M.mut,'patient_idx') && ~isfield(M.mut,'pat_idx'), M.mut=rename_field(M,'patient_idx','pat_idx'); end
demand_fields(M,{'mut','pat'});
if ~isfield(M.pat,'typename') && isfield(M.pat,'ttype_long'), M.pat.typename = M.pat.ttype_long; end
if ~isfield(M.pat,'typename') && isfield(M.pat,'ttype'), M.pat.typename = M.pat.ttype; end
if ~isfield(M.pat,'name') && isfield(M.pat,'indiv'), M.pat.name = M.pat.indiv; end
demand_fields(M.pat,{'name','typename'});
try M.mut = add_simple_fieldnames(M.mut); catch me; end
if ~isfield(M.mut,'newbase'), M.mut.newbase = find_newbase(M.mut); end
demand_fields(M.mut,{'classification','type','newbase','ref_allele'});

fprintf('Preprocessing data for total-rate plot...\n');

% analyze only coding SNPs
M.mut = reorder_struct(M.mut,strcmp(M.mut.classification,'') | strcmp(M.mut.classification,'SNP'));
M.mut = reorder_struct(M.mut,grepi('missense|nonsense|splice|synon|silent|non.?stop|read.?through',M.mut.type,1));
bases = {'A';'C';'G';'T'};
M.mut.from = listmap(M.mut.ref_allele,bases);
M.mut.to = listmap(M.mut.newbase,bases);
M.mut = reorder_struct(M.mut,~isnan(M.mut.from)&~isnan(M.mut.to));

% count mutations per patient
if isfield(M.mut,'patient') && isnumeric(M.mut.patient)
  M.mut.pat_idx = M.mut.patient;
elseif isfield(M.mut,'patient_name')
  M.mut.pat_idx = listmap(M.mut.patient_name,M.pat.name);
elseif iscell(M.mut.patient)
  M.mut.pat_idx = listmap(M.mut.patient,M.pat.name);
else
  error('Don''t know how to map patients!');
end
M.pat.nmut = histc(M.mut.pat_idx,1:slength(M.pat));

% remove patients with zero coding mutations
%M.pat.orig_idx = (1:slength(M.pat))';
%M.pat = reorder_struct(M.pat,M.pat.nmut>0);
%M.mut.pat_idx = mapacross(M.mut.pat_idx,M.pat.orig_idx,(1:slength(M.pat))');

% assume full coverage
M.pat.N = repmat(30e6,slength(M.pat),1);
M.pat.mutrate = M.pat.nmut ./ M.pat.N;
M.pat.lograte = log10(M.pat.mutrate);

% process tumor types
ttype = []; [ttype.name tmp cj] = unique(M.pat.typename);
ttype.median = nan(slength(ttype),1);
ttype.count = nan(slength(ttype),1);
for i=1:slength(ttype)
  idx = find(cj==i); 
  ttype.count(i) = length(idx);
  ttype.median(i) = median(M.pat.lograte(idx));
end
ttype = sort_struct(ttype,'median');
M.pat.ttype_idx = listmap(M.pat.typename,ttype.name);
%if P.include_edits_for_print_version
%  ttype.name = add_newlines_to_long_strings(ttype.name,25);
%else
  ttype.name = regexprep(ttype.name,'lymphocytic','lympho\-cytic');  % (force hyphenation)
  ttype.name = regexprep(ttype.name,'Lung adenocarcinoma','Lung adeno\-carcinoma');
  ttype.name = regexprep(ttype.name,'myelomonocytic','myelo\-monocytic');
  ttype.name = add_newlines_to_long_strings(ttype.name,20);
%end

fprintf('Preprocessing data for spectra boxes...\n');
map = [0 1 3 2,5 0 4 6,6 4 0 5,2 3 1 0]';
if P.include_edits_for_print_version
  classcolors = [1 1 0;0 0.70 0.95;1 0 0;0.1 0.8 0.1;0.7 0.3 0.7;0 0.4 1];
else
  classcolors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0.5 0.3 0.7;0 0.2 0.8];
end
class = {'T->G','T->A','T->C','C->G','C->A','C->T'};
legendparams = {'interpreter','none','fontsize',12};
class = flip(class); map=7-map;
M.mut.frto = ((M.mut.from-1)*4+M.mut.to);
M.mut.class_idx = nansub(map,M.mut.frto);
M.pat.hist = hist2d_fast(M.mut.pat_idx,M.mut.class_idx,1,slength(M.pat),1,length(class));
M.pat.spectrum = bsxfun(@rdivide,M.pat.hist,sum(M.pat.hist,2));

% sort patients within each tumor type

if strcmpi(P.hamburger_sample_sort_method,'rate')

  [tmp ord] = sort(M.pat.mutrate);

elseif strcmp(P.hamburger_sample_sort_method,'spectrum')

  % within each tumor type, cluster by similarity 
  final_ord = cell(slength(ttype),1);
  for tti=1:slength(ttype)
    idx = find(M.pat.ttype_idx == tti);
    q = M.pat.spectrum(idx,:);

    [tmp ord0] = sort(M.pat.mutrate(idx));
    ord1 = dendrogram_get_perm(linkage(pdist(q,'euclidean')),0);
    ord2 = dendrogram_get_perm(linkage(pdist(q,'correlation')),0);
    ord3 = dendrogram_get_perm(linkage(pdist(log10(q),'euclidean')),0);
    ord4 = dendrogram_get_perm(linkage(pdist(log10(q),'correlation')),0);
    [tmp ord_ctr] = sort(sum(M.pat.spectrum(idx,[1]),2));
    [tmp ord_flip] = sort(sum(M.pat.spectrum(idx,[3 5]),2),'descend');

    softmax = -log10(M.pat.hist(idx,:));
    softmax(isinf(softmax)) = nan;
    samplemean = nanmean(softmax,1);
    softmax(isnan(softmax)) = 0;
    softmax = bsxfun(@minus,softmax,samplemean);
    ord5 = dendrogram_get_perm(linkage(pdist(softmax,'euclidean')),0);
    ord6 = dendrogram_get_perm(linkage(pdist(softmax,'correlation')),0);

   % GRAPHICAL COMPARISON OF DIFFERENT ORDERING METHODS
   if 0
    clf;params = {1,'stacked','linestyle','none'};
    subplot(3,3,1); bar(q(ord0,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('mutrate');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    text(0,0.5,[ttype.name{tti} '  '],'fontsize',20,'horizontalalignment','right');set(gcf,'color',[1 1 1]);
    subplot(3,3,2); bar(q(ord_ctr,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('frac(C->T)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    subplot(3,3,3); bar(q(ord_flip,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('frac(C->G + A->T)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    subplot(3,3,4); bar(q(ord1,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('euclidean(frac)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    subplot(3,3,7); bar(q(ord2,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('correlation(frac)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    subplot(3,3,5); bar(q(ord3,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('euclidean(logfrac)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    subplot(3,3,8); bar(q(ord4,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('correlation(logfrac)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    subplot(3,3,6); bar(q(ord5,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('euclidean(softmax)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
    subplot(3,3,9); bar(q(ord6,:),params{:}); xlim([0.5 length(idx)+0.5]); ylim([0 1]);title('correlation(softmax)');
    set(gca,'ydir','reverse','xtick',[],'ytick',[]);box on; colormap(classcolors);
   end
    final_ord{tti} = idx(ord_flip);
  end    
  ord = cat(1,final_ord{:});
end

M.pat.orig_idx = (1:slength(M.pat))';
M.pat = reorder_struct(M.pat,ord);
%M.mut.pat_idx = mapacross(M.mut.pat_idx,M.pat.orig_idx,(1:slength(M.pat))');
% (don't need to do this, because we never consult M.mut again)

% DRAW FIGURE

figure(1); clf; set(gcf,'position',[30 149 1402 665]);
set(gcf,'paperpositionmode','auto','papersize',[11 8.5])

% TOTAL RATE PLOT
left = 0.1; width = 0.87; height = 0.52; ydivide = 0.43;
subplot('position',[left ydivide width height]);
if P.include_edits_for_print_version
  P = impose_default_value(P,'fontsize',15);
end
ph = catplot(M.pat.lograte,M.pat.ttype_idx,zeros(slength(M.pat),3),ttype.name,0.5,'o',NaN,[-8 -3],P);
set(gca,'tickdir','out','ticklength',[0.003 0.002]);
set(gca,'ytick',-8:-3,'yticklabel',{'0.01','0.1','1','10','100','1000'});
ylabel('Somatic mutation frequency (/Mb)','fontsize',20);
set(gcf,'color',[1 1 1]);set(cat(2,ph{:}),'MarkerSize',3);
for i=1:slength(ttype)
  line(i+[-0.3 0.3],ttype.median(i)*[1 1],'color',[0.5 0.5 0.5]);
  ct = num2str(ttype.count(i)); if i==1, ct = ['n=' ct]; end
  text(i,-2.8,ct,'fontsize',12,'horizontalalignment','center');
end

% SPECTRA BOXES
eachwidth = width / slength(ttype);
specheight = 0.12;
for i=1:slength(ttype)
  subplot('position',[left + (i-1)*(eachwidth), 0.03, eachwidth*0.9, specheight]);
  idx=find(M.pat.ttype_idx==i);
  if length(idx)>1
    bar(M.pat.spectrum(idx,:),1,'stacked');
  elseif length(idx)==1  % workaround to force correct display of single column
    q = M.pat.spectrum(idx,:); q(2,:) = 0;
    bar(q,1,'stacked');
  end
  shading flat;ylim([0 1]);xlim([0.5 length(idx)+0.5]);
  set(gca,'YTick',[],'XTick',[],'ydir','reverse');box on
  line([0.5 length(idx)+0.5 length(idx)+0.5 0 0],[0 0 1 1 0],'Color',[0 0 0]);
end
set(gcf,'color',[1 1 1]);colormap(classcolors);
legend(class,'position',[0.02 0.05 0.06 0.15],legendparams{:}); legend boxoff

if ~isempty(P.hamburger_plot_output_filename)
  print_to_file(P.hamburger_plot_output_filename);
end


