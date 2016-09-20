function R = analyze_gmzc65_data(Nn,outname)

Z = load_struct('/xchip/cga1/lawrence/db/gmzc65/categs.txt');
tmp = parse(Z.name,'(.*):(.*):(.*):(.*)$',{'goodbad','expr','zone','categ'});
Z = merge_structs({Z,tmp});

% remove "bad" regions and Ns
idx =grepv('bad|any N',Z.name,1);
Z = reorder_struct(Z,idx);
Nn = Nn(idx,:);

% parse 'categ'
tmp = parse(Z.categ,'^([ACGT]) in (.*)$',{'base','context'});
Z = merge_structs({Z,tmp});

% collapse by strand
idx = grep('G|T',Z.base,1);
Z.base(idx) = rc(Z.base(idx));
Z.context(idx) = rc(Z.context(idx));
z = Nn(idx,5:-1:2); Nn(idx,2:5) = z;
tmp = stringsplice([Z.expr Z.zone Z.base Z.context],1,':');
[u ui uj] = unique(tmp,'first');
out = zeros(length(u),5);
for i=1:length(uj), out(uj(i),:)=out(uj(i),:)+Nn(i,:); end
Z = reorder_struct(Z,ui);
Nn = out;

% calculate rates

q = {'unknown','zero_expression','nonzero_expression','expression_5pct',...
   'expression_25pct','expression_50pct','expression_75pct'};
qlabel = q;
z = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
zlabel = {'IGR','intron','UTR/pro','exon','total'};
clabel = {'CpG transition','CpG transversion', 'other C:G transition', 'other C:G transversion', 'A:T transition', ...
          'A:T transversion', 'total'};
Z.ctype = 3*ones(slength(Z),1);
idx = grep('C',Z.base,1); Z.ctype(idx) = 2;
idx = idx(grep('_G',Z.context(idx),1)); Z.ctype(idx) = 1;
R = []; R.categ=[]; R.zone=[]; R.context=[]; R.Ncov=[]; R.n_muts=[]; R.mutrate=[]; r=1;
R.stdev=[]; R.ci_low=[]; R.ci_high=[];
for qi=1:length(qlabel), idx1 = find(strcmp(q{qi},Z.expr));
  for ci=1:length(clabel)
    if ~strcmp(clabel{ci},'total')
      idx2 = idx1(Z.ctype(idx1)==ceil((ci-0.1)/2));
      if ci==1 || ci==3, cols = 5; elseif ci==2 || ci==4, cols = 2:4; else cols = 2:5; end
      if ci==5, cols = 4; elseif ci==6, cols = [2 3 5]; end
    else
      idx2 = idx1; cols=2:5;
    end
    for zi=1:length(zlabel), idx3 = idx2(grep(z{zi},Z.zone(idx2),1));
      R.categ{r,1} = qlabel{qi};  R.zone{r,1} = zlabel{zi};  R.context{r,1} = clabel{ci};
      R.Ncov(r,1) = sum(Nn(idx3,1));  R.n_muts(r,1) = fullsum(Nn(idx3,cols));
      [rate confint] = binofit(R.n_muts(r,1),R.Ncov(r,1));
      R.mutrate(r,1) = 1e6*rate;
      R.ci_low(r,1) = 1e6*confint(1);
      R.ci_high(r,1) = 1e6*confint(2);      
      R.stdev(r,1) = (R.ci_high(r,1)-R.mutrate(r,1))/1.96;
      r=r+1;
end,end,end            

save_struct(R,outname);
