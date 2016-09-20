function R = analyze_globcov_data(Nn,categdir,outstem,P)
% R = analyze_globcov_data(Nn,categdir,outstem,P)
%
% input:
%        Nn = N n->A n->C n->G n->T
%
%
% general method: allows user to arbitrarily configure:
%
%     x  categories to exclude up-front (e.g. "bad", "any N")
%
%     q  (major division, shown as columns)
%     z  (minor division, shown as rows within each base-context category)
%     c  (base-context categories)
%         currently:  only two options
%                     1  (CpG/CG/AT)x(transition/transversion)
%                     2  melanoma categories
%         eventually: will accept a K struct from automatic mutation discovery
%
% Mike Lawrence 2010-06-16

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'exclude_grep',[]);
P=impose_default_value(P,'major_grep','*required*');
P=impose_default_value(P,'major_label','*required*');
P=impose_default_value(P,'minor_grep','*required*');
P=impose_default_value(P,'minor_label','*required*');
P=impose_default_value(P,'output_format',2);
P=impose_default_value(P,'use_melanoma_categories',false);

if categdir(1)=='/'
  fname = [categdir '/categs.txt'];
else
  fname = ['/xchip/cga1/lawrence/db/' categdir '/categs.txt'];
end
Z = load_struct(fname);
if size(Nn,1)~=slength(Z), error('size(Nn,1)~=slength(Z)'); end
if size(Nn,2)~=5, error('size(Nn,2)~=5'); end
Z.Nn = Nn; clear Nn;

% remove "bad" regions, noncoding regions, non-exons, and Ns
if ~isempty(P.exclude_grep)
  fprintf('Excluding: %s\n',P.exclude_grep);
  idx = grepv(P.exclude_grep,Z.name,1);
  Z = reorder_struct(Z,idx);
end

% find base and context
fprintf('Finding base and context\n');
tmp = parse(Z.name,'([ACGT]) in ([^:]*)',{'base','context'});
if any(cellfun('isempty',tmp.context) | cellfun('isempty',tmp.base)), error('Error determining base and context'); end
Z = merge_structs({Z,tmp});

% for base=G or T, reverse-complement the base, context, and mutation counts
fprintf('Collapsing strands\n');
idx = grep('G|T',Z.base,1);
Z.base(idx) = rc(Z.base(idx));
Z.context(idx) = rc(Z.context(idx));
tmp = Z.Nn(idx,5:-1:2); Z.Nn(idx,2:5) = tmp;

%keyboard
% identity territories of the base-context categories
if P.use_melanoma_categories
  clabel = {...
    'py-surrounded_C_transition','py-surrounded_C_transversion',...
    'py-adjacent_C_transition','py-adjacent_C_transversion',...
    'non-py-adjacent_CpG_transition','non-py-adjacent_CpG_transversion',...
    'other_C_transition', 'other_C_transversion',...
    'A_transition','A_transversion',...
    'total'};
  Z.ctype = 5*ones(slength(Z),1);
  idx = grep('C',Z.base,1); Z.ctype(idx) = 4;
  idx2 = idx(grep('_G',Z.context(idx),1)); Z.ctype(idx2) = 3;
  idx2 = idx(grep('[CT]_|_[CT]',Z.context(idx),1)); Z.ctype(idx2) = 2;
  idx2 = idx(grep('[CT]_[CT]',Z.context(idx),1)); Z.ctype(idx2) = 1;
elseif P.use_esophageal_categories
  clabel = {...
      'CpG_to_A_or_T', 'ApG_to_C_or_G', ...
      'Np(A/C/T)-nonflip', 'flip', 'total'}
  Z.ctype = 4*ones(slength(Z),1);
  idx = grep('_A|_C|_T', Z.base, 1);
  Z.ctype(idx) = 3; 
  idx = grep('C', Z.base, 1); 
  idx = idx(grep('_G',Z.context(idx),1));
  Z.ctype(idx) = 1;
  idx = grep('A', Z.base, 1);
  idx = idx(grep('_G',Z.context(idx),1));
  Z.ctype(idx) = 2;
else
  clabel = {...
     'CpG_transition','CpG_transversion',...
     'other_C_transition','other_C_transversion',...
     'A_transition','A_transversion',...
     'total'};
  Z.ctype = 3*ones(slength(Z),1);
  idx = grep('C',Z.base,1); 
  Z.ctype(idx) = 2;
  idx = idx(grep('_G',Z.context(idx),1)); 
  Z.ctype(idx) = 1;
end

% compute rates
R = []; r=1;
fprintf('Computing rates\n');
for maji=1:length(P.major_label), idx1 = grep(P.major_grep{maji},Z.name,1);
  for ci=1:length(clabel)
    if strcmp(clabel{ci},'total')
      idx2 = idx1; cols=2:5;
    else
      idx2 = idx1(Z.ctype(idx1)==ceil((ci-0.1)/2));
      if P.use_melanoma_categories
        if ci==1||ci==3||ci==5||ci==7, cols=5; elseif ci==2||ci==4||ci==6||ci==8, cols=2:4; else cols=2:5; end
        if ci==9, cols=4; elseif ci==10; cols=[2 3 5]; end
      elseif P.use_esophageal_categories 
        if ci == 1, cols = [2, 5]; elseif ci == 2, cols = [3, 4]; elseif ci == 3, cols = 2:5, else cols = 2:5;, end
      else
        if ci==1||ci==3, cols = 5; elseif ci==2||ci==4, cols = 2:4; else cols = 2:5; end
        if ci==5, cols=4; elseif ci==6, cols=[2 3 5]; end
      end
    end
    for mini=1:length(P.minor_label), 
      idx3 = idx2(grep(P.minor_grep{mini},Z.name(idx2),1));
      major_field{r,1} = P.major_label{maji};  minor_field{r,1} = P.minor_label{mini};  
      R.base_context{r,1} = clabel{ci};
      R.Ncov(r,1) = sum(Z.Nn(idx3,1));  
      R.n_muts(r,1) = fullsum(Z.Nn(idx3,cols));
      [rate confint] = binofit(R.n_muts(r),R.Ncov(r));
      R.mutrate(r,1) = 1e6*rate;
      R.ci_low(r,1) = 1e6*confint(1);
      R.ci_high(r,1) = 1e6*confint(2);      
      R.stdev(r,1) = (R.ci_high(r)-R.mutrate(r))/1.96;
      r=r+1;
end,end,end            
R = setfield(R,P.major_title,major_field);
R = setfield(R,P.minor_title,minor_field);
R = order_fields_first(R,{P.major_title,'base_context',P.minor_title,'Ncov','n_muts','mutrate','stdev','ci_low','ci_high'});

% write output
fprintf('Writing output\n');
outname = [outstem '.analysis.mat'];
save(outname,'R');
outname = [outstem '.analysis.txt'];

if P.output_format==1
  save_struct(R,outname);

elseif P.output_format==2
  out = fopen(outname,'wt');
  flds = fieldnames(R); 
  majorfld = P.major_title; minorflds = {'base_context',P.minor_title};
  otherflds = flds(~ismember(flds,union(majorfld,minorflds)));
  otherflds_formats = {'%d','%d','%0.3f','%0.3f','%0.3f','%0.3f'};
  tmp = getfield(R,majorfld); [u ui uj] = unique(tmp);  majorfld_members = tmp(sort(ui));
  majf = getfield(R,majorfld);
  minf = cell(length(minorflds),1); for i=1:length(minorflds), minf{i} = getfield(R,minorflds{i}); end
  othf = cell(length(otherflds),1); for i=1:length(otherflds), othf{i} = getfield(R,otherflds{i}); end
  for i=1:length(minorflds), fprintf(out,'\t'); end
  for i=1:length(majorfld_members)
    fprintf(out,'\t');
    for j=1:length(otherflds)
      if j==round(length(otherflds)/2), fprintf(out,'%s',majorfld_members{i}); end
      fprintf(out,'\t');
    end
  end
  fprintf(out,'\n');
  for i=1:length(minorflds), fprintf(out,'%s\t',minorflds{i}); end
  for i=1:length(majorfld_members)
    fprintf(out,'\t');
    for j=1:length(otherflds), fprintf(out,'%s\t', otherflds{j}); end
  end
  fprintf(out,'\n');
  for row=1:slength(R)
    if ~strcmp(majf{row},majorfld_members{1}), break; end
    if row==1 || ~strcmp(minf{1}{row},minf{1}{row-1}), fprintf(out,'\n'); end
    for i=1:length(minorflds), fprintf(out,'%s\t',minf{i}{row}); end
    for i=1:length(majorfld_members)
      idx = find(strcmp(majf,majorfld_members{i}));
      for j=1:length(minorflds), idx=idx(strcmp(minf{j}(idx),minf{j}{row})); end
      fprintf(out,'\t');
      for j=1:length(otherflds), fprintf(out,[otherflds_formats{j} '\t'], othf{j}(idx)); end
    end
    fprintf(out,'\n');
  end
  fclose(out);
else
  error('Unknown P.output_format');
end

