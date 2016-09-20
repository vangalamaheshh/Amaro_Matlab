function R = analyze_gcszp_data(Nn,outname,P)

if ~exist('P','var'),P=[];end
P=impose_default_value(P,'use_melanoma_categories',false);
P=impose_default_value(P,'output_format',1);

Z = load_struct('/xchip/cga1/lawrence/db/gcsz29p/categs.txt');  % (same as gcsz3p)
tmp = parse(Z.name,'(.*):([ACGT]) in (.*)$',{'categ','base','context'});
Z = merge_structs({Z,tmp});

% remove "bad" regions and Ns
idx =grepv('^(bad|any N)',Z.name,1);
Z = reorder_struct(Z,idx);
Nn = Nn(idx,:);

% collapse by strand
idx = grep('G|T',Z.base,1);
Z.base(idx) = rc(Z.base(idx));
Z.context(idx) = rc(Z.context(idx));
z = Nn(idx,5:-1:2); Nn(idx,2:5) = z;
tmp = stringsplice([Z.categ Z.base Z.context],1,':');
[u ui uj] = unique(tmp,'first');
out = zeros(length(u),5);
for i=1:length(uj), out(uj(i),:)=out(uj(i),:)+Nn(i,:); end
Z = reorder_struct(Z,ui);
Nn = out;

% compute rates

q = {':cons','noncons','cons'}; qlabel = {'cons','noncons','total'};
z = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
zlabel = {'IGR','intron','UTR/pro','exon','total'};

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
else
  clabel = {...
     'CpG_transition','CpG_transversion',...
     'other_C_transition','other_C_transversion',...
     'A_transition','A_transversion',...
     'total'};
  Z.ctype = 3*ones(slength(Z),1);
  idx = grep('C',Z.base,1); Z.ctype(idx) = 2;
  idx = idx(grep('_G',Z.context(idx),1)); Z.ctype(idx) = 1;
end

R = []; R.categ=[]; R.context=[]; R.zone=[]; R.Ncov=[]; R.n_muts=[]; R.mutrate=[]; r=1;
R.stdev=[]; R.ci_low=[]; R.ci_high=[];
for qi=1:length(qlabel), idx1 = grep(q{qi},Z.categ,1);
  for ci=1:length(clabel)
    if strcmp(clabel{ci},'total')
      idx2 = idx1; cols=2:5;
    else
      idx2 = idx1(Z.ctype(idx1)==ceil((ci-0.1)/2));
      if P.use_melanoma_categories
        if ci==1||ci==3||ci==5||ci==7, cols=5; elseif ci==2||ci==4||ci==6||ci==8, cols=2:4; else cols=2:5; end
        if ci==9, cols=4; elseif ci==10; cols=[2 3 5]; end
      else
        if ci==1||ci==3, cols = 5; elseif ci==2||ci==4, cols = 2:4; else cols = 2:5; end
        if ci==5, cols=4; elseif ci==6, cols=[2 3 5]; end
      end
    end
    for zi=1:length(zlabel), idx3 = idx2(grep(z{zi},Z.categ(idx2),1));
      R.categ{r,1} = qlabel{qi};  R.zone{r,1} = zlabel{zi};  R.context{r,1} = clabel{ci};
      R.Ncov(r,1) = sum(Nn(idx3,1));  R.n_muts(r,1) = fullsum(Nn(idx3,cols));
      [rate confint] = binofit(R.n_muts(r,1),R.Ncov(r,1));
      R.mutrate(r,1) = 1e6*rate;
      R.ci_low(r,1) = 1e6*confint(1);
      R.ci_high(r,1) = 1e6*confint(2);      
      R.stdev(r,1) = (R.ci_high(r,1)-R.mutrate(r,1))/1.96;
      r=r+1;
end,end,end            

if P.output_format==1
  save_struct(R,outname);

elseif P.output_format==2
  out = fopen(outname,'wt');
  flds = fieldnames(R); 
  majorfld = 'categ'; minorflds = {'context','zone'};
  otherflds = flds(~ismember(flds,union(majorfld,minorflds)));
  otherflds_formats = {'%d','%d','%0.3f','%0.3f','%0.3f','%0.3f'};
  tmp = getfield(R,majorfld); [u ui uj] = unique(tmp);  majorfld_members = tmp(ui);
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

