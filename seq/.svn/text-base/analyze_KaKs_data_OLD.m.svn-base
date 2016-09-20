function [R H] = analyze_KaKs_data(Ntot,mutrates,outstem)
% Ntot = number of covered bases in each of the 1690 "gc65e" categories
% mutrates = [CpG transition, CpG transversion, other C:G transition, other C:G transversion, A:T transition, A:T transversion]
%
% returns:
%  R = distribution of expected nonsilent/silent ratio
%  H = distributions of expected nonsilent and silent mutation counts

Ntot = as_column(Ntot);
mutrates = as_column(mutrates);

if size(Ntot) ~= [1690 1], error('invalid Ntot'); end
if size(mutrates) ~= [6 1], error('invalid mutrates'); end
if any(mutrates<1e-10 | mutrates>1e-4), error('mutrates out of range'); end

Z = load_struct('/xchip/cga1/lawrence/db/gc65e/categs.txt');
Z.Ntot = Ntot;

% remove "bad" regions, noncoding regions, and Ns
idx = grepv('any N|bad|noncoding',Z.name,1);
Z = reorder_struct(Z,idx);

% parse names
tmp = parse(Z.name,'^(good|bad):([ACGT]) in ([ACGT])_([ACGT]):(.*)$',...
  {'goodbad','base','left','right','effect'});
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
newZ.Ntot = zeros(slength(newZ),1);
for i=1:length(uj), newZ.Ntot(uj(i)) = newZ.Ntot(uj(i)) + Z.Ntot(i); end
Z=newZ;

% Remaining steps:
% (1) make table X:  N newbase silent/nonsilent rate_to_newbase (omit rows where N=0)
% (2) separate into two tables Xs and Xn (silent + nonsilent)
% (3) convolute silent table to get distribution of expected number of silent mutations
% (4) convolute nonsilent table to get distribution of expected number of nonsilent mutations
% (5) divide the two distributions to get distribution of expected nonsilent/silent ratio

% expand table into triplicate
base = 'ACGT';
Y = cell(4,1);
for i=1:4
  Y{i} = Z;
  Y{i}.newbase = repmat({base(i)},slength(Y{i}),1);
end
Y = concat_structs(Y);
Y = reorder_struct(Y,~strcmp(Y.base,Y.newbase));

% omit rows where Ntot==0
Y = reorder_struct(Y,Y.Ntot>0);

% assign rates
transition(base) = 'GTAC';
Y.istransition = false(slength(Y),1);
for i=1:slength(Y), Y.istransition(i) = (Y.newbase{i}==transition(Y.base{i})); end
Y.rate = mutrates(2*Y.ctype - Y.istransition);
Y.rate(~Y.istransition) = Y.rate(~Y.istransition)/2;   % (transversion rates need to be halved)

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

X = cell(2,1);   % (nonsilent,silent)
X{1} = reorder_struct(Y,~Y.silent);
X{2} = reorder_struct(Y,Y.silent);

% compute expected number of mutations for each
nexp = nan(2,1); for i=1:2, nexp(i) = round(sum(X{i}.Ntot .* X{i}.rate)); end
% compute maximum number of mutations to consider
maxmuts = nexp * 2;
tailtrim = 1e-30;

% convolute each table
type = {'nonsilent','silent'};
for i=1:2
  fprintf('Computing histograms for %s\n',type{i});
  H{i} = cell(slength(X{i}),1);
  for j=1:slength(X{i})
    H{i}{j} = binopdf(0:maxmuts,X{i}.Ntot(j),X{i}.rate(j));
  end
  fprintf('Convoluting histograms for %s\n',type{i});
  H{i} = batch_convolute(H{i});
  % trim tail
  trimpoint = find(H{i}>tailtrim,1,'last');
  H{i}=H{i}(1:trimpoint);
end

% divide the two distributions to get the distribution of the ratio
fprintf('Computing distribution of expected nonsilent/silent ratio\n');
maxratio = 10;
binsize = 0.05;
numbins = ceil(maxratio/binsize)+1;
R = zeros(numbins+1,1);   % add one more bin at end for "overflow"
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

[mx idx] = max(R);
mle = (idx-1)*binsize;

out = fopen([outstem '.MLE.txt'],'wt');
fprintf('Maximum-likelihood estimates:\n');
fprintf('\tNumber of nonsilent mutations:   %d\n', nexp(1));
fprintf('\tNumber of silent mutations:      %d\n', nexp(2));
fprintf('\tKa/Ks = nonsilent/silent ratio:  %f\n', mle);

fprintf(out,'Maximum-likelihood estimates:\n');
fprintf(out,'\tNumber of nonsilent mutations:   %d\n', nexp(1));
fprintf(out,'\tNumber of silent mutations:      %d\n', nexp(2));
fprintf(out,'\tKa/Ks = nonsilent/silent ratio:  %f\n', mle);

fclose(out);




keyboard
