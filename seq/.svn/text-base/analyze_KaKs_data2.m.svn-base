function [R H rbins pval] = analyze_KaKs_data2(Nn,categdir,outstem,P)
% [R H rbins] = analyze_KaKs_data2(Nn,categdir,outstem)
% given:
%  Nn = rows: categories as in categdir; columns: 1=N 2=n->A 3=n->C 4=n->G 5=n->T
%  categdir = should include effect and "A in C_G" context
% returns:
%  R = distribution of expected nonsilent/silent ratio
%  H = distributions of expected nonsilent and silent mutation counts
%  rbins = bin centers of each position in R

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'title',[]);
P=impose_default_value(P,'method',1);
P=impose_default_value(P,'remove_bad',true);

if categdir(1)=='/'
  fname = [categdir '/categs.txt'];
else
  fname = ['/xchip/cga1/lawrence/db/' categdir '/categs.txt'];
end
Z = load_struct(fname);
if size(Nn,1)~=slength(Z), error('size(Nn,1)~=slength(Z)'); end
Z.Nn = Nn; clear Nn;

% remove noncoding regions, non-exons, and Ns
fprintf('Removing noncoding regions, non-exons, and Ns\n');
idx = grepv('any N|noncoding|intron|IGR|UTR',Z.name,1);
Z = reorder_struct(Z,idx);
% remove "bad" regions
if P.remove_bad
  fprintf('Removing bad regions\n');
  idx = grepv('bad',Z.name,1);
  Z = reorder_struct(Z,idx);
end

% parse names
tmp = parse(Z.name,'^(good|bad):(.*:)?([ACGT]) in ([ACGT])_([ACGT]):(.*)$',...
  {'goodbad','other','base','left','right','effect'});
Z = merge_structs({Z,tmp});

% collapse 64 categories to CpG, other C:G, A:T
Z.ctype = 3*ones(slength(Z),1);
Z.ctype(grep('C|G',Z.base,1)) = 2;
Z.ctype((strcmp(Z.base,'C')&strcmp(Z.right,'G'))|...
  strcmp(Z.base,'G')&strcmp(Z.left,'C')) = 1;
context = {'CpG';'other C:G';'A:T'};
Z.context = context(Z.ctype);
Z.name = stringsplice([Z.base Z.context Z.effect],1,':');
Z = rmfield(Z,{'right','left'});
[u ui uj] = unique(Z.name);
newZ = reorder_struct(Z,ui);
newZ.Nn = zeros(size(newZ.Nn));
for i=1:length(uj)
  newZ.Nn(uj(i),:) = newZ.Nn(uj(i),:) + Z.Nn(i,:);
end
Z=newZ;

base = 'ACGT'; transition(base) = 'GTAC';

% expand table into triplicate
Y = cell(4,1);
for i=1:4
  Y{i} = Z;
  Y{i}.newbase = repmat({base(i)},slength(Y{i}),1);
  Y{i}.N = Y{i}.Nn(:,1);
  Y{i}.n = Y{i}.Nn(:,1+i);
end
Y = concat_structs(Y);
Y = reorder_struct(Y,~strcmp(Y.base,Y.newbase));
Y = rmfield(Y,'Nn');
Y.istransition = false(slength(Y),1);
for i=1:slength(Y), Y.istransition(i) = (Y.newbase{i}==transition(Y.base{i})); end

% assign whether each row is silent or nonsilent
Y.silent = false(slength(Y),1);
[u ui uj] = unique(Y.effect);
for i=1:length(u)
  idx = find(uj==i);
  if strcmp(u{i},'any change is nonsilent')
    continue;
  elseif strcmp(u{i},'any change is silent')
    Y.silent(idx)=true;
  elseif contains(u{i},'/')   % two possible silent changes
    tmp = parse(u{i},'change to ([ACGT])/([ACGT]) is silent',{'s1','s2'});
    if isempty(tmp.s1{1}) || isempty(tmp.s2{1}), error('Invalid effect: %s',u{i}); end
    idx = idx(grep([tmp.s1{1} '|' tmp.s2{1}],Y.newbase(idx),1));
    Y.silent(idx)=true;
  else   % one possible silent change
    tmp = parse(u{i},'change to ([ACGT]) is silent',{'s'});
    if isempty(tmp.s{1}), error('Invalid effect: %s',u{i}); end
    idx = idx(grep(tmp.s{1},Y.newbase(idx),1));
    Y.silent(idx)=true;
  end
end

nobs_nonsilent = sum(Y.n(~Y.silent));
nobs_silent = sum(Y.n(Y.silent));
obs_ratio = nobs_nonsilent / nobs_silent;

if P.method==1
  % first method developed:
  %
  % given tabulated coverage and mutation counts per context category,
  % (1) computes mutation rates per context category
  % (2) computes distribution of expected number of silent mutations
  % (3) computes distribution of expected number of nonsilent mutations
  % (4) divides (3)/(2) to get distribution of expected NS/S ratio

  % compute mutation rates from the data:
  %   (CpG transit, CpG transver, other CG transit, other CG transver, AT transit, AT transver)
  
  z = zeros(length(context),1); Ncov=z; ntransit=z; ntransver=z; 
  for i=1:length(context)
    idx = find(strcmp(Z.context,context{i}));
    for j=1:length(idx), k=idx(j);
      Ncov(i)=Ncov(i)+Z.Nn(k,1);
      for b=1:4
        if strcmp(base(b),Z.base{k}), continue;  % non-mutation
        elseif strcmp(base(b),transition(Z.base{k})), ntransit(i)=ntransit(i) + Z.Nn(k,b+1);  % transition
        else ntransver(i)=ntransver(i) + Z.Nn(k,b+1);
  end,end,end,end
  mutrates = nan(length(context)*2,1);
  for i=1:length(context), mutrates(i*2+[-1 0]) = [ntransit(i) ntransver(i)] / Ncov(i); end
  
  % Remaining steps:
  % (1) make table X:  N newbase silent/nonsilent rate_to_newbase (omit rows where N=0)
  % (2) separate into two tables Xs and Xn (silent + nonsilent)
  % (3) convolute silent table to get distribution of expected number of silent mutations
  % (4) convolute nonsilent table to get distribution of expected number of nonsilent mutations
  % (5) divide the two distributions to get distribution of expected nonsilent/silent ratio
 
  % omit rows where Ncov==0
  Y = reorder_struct(Y,Y.N>0);
  
  % assign rates
  Y.rate = mutrates(2*Y.ctype - Y.istransition);
  Y.rate(~Y.istransition) = Y.rate(~Y.istransition)/2;   % (transversion rates need to be halved)
  
  % separate into nonsilent and silent
  X = cell(2,1);   % (nonsilent,silent)
  X{1} = reorder_struct(Y,~Y.silent);
  X{2} = reorder_struct(Y,Y.silent);
  
  % compute expected number of mutations for each
  nexp = nan(2,1); for i=1:2, nexp(i) = round(sum(X{i}.N .* X{i}.rate)); end
  % compute maximum number of mutations to consider
  maxmuts = nexp * 2;
  tailtrim = 1e-30;
  
  % convolute each table
  type = {'nonsilent','silent'};
  for i=1:2
    fprintf('Computing histograms for %s\n',type{i});
    H{i} = cell(slength(X{i}),1);
    for j=1:slength(X{i})
      H{i}{j} = binopdf(0:maxmuts,X{i}.N(j),X{i}.rate(j));
    end
    fprintf('Convoluting histograms for %s\n',type{i});
    H{i} = batch_convolute(H{i});
    % trim tail
    trimpoint = find(H{i}>tailtrim,1,'last');
    H{i}=H{i}(1:trimpoint);
  end
  
  % divide the two distributions to get the distribution of the ratio
  fprintf('Computing distribution of expected nonsilent/silent ratio\n');
  maxratio = 50;
  if sum(nexp)>500, binsize = 0.05;
  elseif sum(nexp)>100, binsize = 0.2;
  else binsize = 0.4;
  end
  
  numbins = ceil(maxratio/binsize)+1;
  R = zeros(numbins+1,1);   % add one more bin at end for "overflow"
  rbins = [(1:numbins)'*binsize;inf];
  for i=1:length(H{1}), if ~mod(i,100), fprintf('%d/%d ',i,length(H{1})); end
    for j=1:length(H{2})
      p = H{1}(i) * H{2}(j);
      n_nonsilent = i-1;
      n_silent = j-1;
      r = n_nonsilent/n_silent;
      if isnan(r) || isinf(r) || r>maxratio
        binno = numbins+1;
      else
        binno = round(r/binsize)+1;
      end
      R(binno) = R(binno) + p;
  end,end, fprintf('\n');

  nexp_nonsilent = nexp(1);
  nexp_silent = nexp(2);
  
elseif P.method==2
  % second method developed:
  %
  % given tabulated coverage and mutation counts per context category x effect category
  % (1) in each context category, distributes mutation counts between the effect categories
  % (2) convolutes distributions to build distribution of NS count
  % (3) trivially derives distribution of S count, and NS/S ratio
  %
  % Note: preserves mutation context category, base, and newbase

  % Make new table W: context base newbase n_obs N_where_it_would_be_silent N_where_it_would_be_nonsilent

  [c ci cj] = unique(Y.context);
  [b bi bj] = unique(Y.base);
  [n ni nj] = unique(Y.newbase);
  W = []; widx = 1;
  for cidx=1:length(c), for bidx=1:length(b), for nidx=1:length(n)
    idx = find(cj==cidx & bj==bidx & nj==nidx);
    W.context{widx,1} = c{cidx};
    W.base{widx,1} = b{bidx};
    W.newbase{widx,1} = n{nidx};
    W.n(widx,1) = sum(Y.n(idx));
    W.Nsil(widx,1) = sum(Y.N(idx(Y.silent(idx))));
    W.Nnon(widx,1) = sum(Y.N(idx(~Y.silent(idx))));
    widx=widx+1;
  end,end,end

  % remove entries with no observed mutations
  W = reorder_struct(W,W.n>0);
  W.fracNnon = W.Nnon./(W.Nnon+W.Nsil);

  % calc NS = distribution of number of nonsilent mutations, for each entry
  nw = slength(W);
  NS = cell(nw,1);   % NS{i}(1)=0 nonsilent NS{i}(2) = 1 nonsilent ...  NS{i}{n+1} = n nonsilent
  for i=1:nw
    NS{i} = binopdf(0:W.n(i),W.n(i),W.fracNnon(i));
  end
    
  % convolute the distributions
  fprintf('Convoluting distributions\n');
  NS = batch_convolute(NS);
  
  % make H (nonsilent, silent)
  S = flipud(NS);
  H = {NS, S};

  [mx idx] = max(NS);
  nexp_nonsilent = idx-1;
  
  [mx idx] = max(S);
  nexp_silent = idx-1;
  
  % make R
  R = NS;
  rbins = (0:sum(W.n))./(sum(W.n):-1:0);

else
  error('uknown P.method');
end

% find MLE for ratio
[mx midx] = max(R);
exp_ratio = rbins(midx);

% find stdev for ratio
sm = nan(midx,1);
st = nan(midx,1);
for i=0:midx-1
  if midx+i>length(R), break; end
  sm(i+1) = sum(R(midx-i:midx+i));
  st(i+1) = (rbins(midx+i)-rbins(midx-i))/2;
end
[df idx] = min(abs(sm-0.682));
if (df<=0.2), stdev_exp_ratio = st(idx); else stdev_exp_ratio = inf; end

% p-value for observed NS/S ratio (1-sided, for having higher ratio)
pval = sum(R(rbins>=obs_ratio));

out = fopen([outstem '.MLE.txt'],'wt');

ffprintf(out,'Maximum-likelihood estimates:\n');
ffprintf(out,'\tNumber of nonsilent mutations:   %d\n', nexp_nonsilent);
ffprintf(out,'\tNumber of silent mutations:      %d\n', nexp_silent);
ffprintf(out,'\tKa/Ks = nonsilent/silent ratio:  %f +- %f (stdev)\n', exp_ratio,stdev_exp_ratio);

ffprintf(out,'\nObserved values:\n');
ffprintf(out,'\tNumber of nonsilent mutations:   %d\n', nobs_nonsilent);
ffprintf(out,'\tNumber of silent mutations:      %d\n', nobs_silent);
ffprintf(out,'\tKa/Ks = nonsilent/silent ratio:  %f\n', obs_ratio);
ffprintf(out,'\tp-value (1-sided, excess NS):    %f\n', pval);

fclose(out);

figure(1);clf

subplot(3,1,1); plot(H{1});
[a b] = max(H{1}); xlim([1 min(b*2,length(H{1}))]);
xlabel('number of nonsilent mutations');
ylabel('probability'); set(gca,'ytick',[]);
if ~isempty(P.title), title(P.title,'interpreter','none'); end
line(nobs_nonsilent*[1 1],ylim,'color',[1 0 0]);

subplot(3,1,2); plot(H{2});
[a b] = max(H{2}); xlim([1 min(b*2,length(H{2}))]);
xlabel('number of silent mutations');
ylabel('probability'); set(gca,'ytick',[]);
line(nobs_silent*[1 1],ylim,'color',[1 0 0]);

subplot(3,1,3); plot(rbins,R);
[a b] = max(R); xlim([1 min(rbins(b)*2,length(R))]);
xlabel('nonsilent/silent ratio');
ylabel('probability'); set(gca,'ytick',[]);
line(obs_ratio*[1 1],ylim,'color',[1 0 0]);

print_to_file([outstem '.png']);

