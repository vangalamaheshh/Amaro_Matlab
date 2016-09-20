function x = load_and_trim_454_data(fname,P)
% Mike Lawrence 2010-06-08

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'leaderseq',[]);
P = impose_default_value(P,'trailerseq',[]);
P = impose_default_value(P,'max_extra_trim',20);
P = impose_default_value(P,'trim_extra_leaders',false);
P = impose_default_value(P,'trailerseq_len_required',30);
P = impose_default_value(P,'trailerseq_score_required',40);
P = impose_default_value(P,'min_readlength_after_trim',200);
P = impose_default_value(P,'soft_trim',false);

if iscell(fname)
  % multi-region mode
  nr = length(fname);
  x = cell(nr,1);
  for i=1:nr
    fprintf('Loading %s\n', fname{i});
    x{i} = load_and_trim_454_data(fname{i},P);
    x{i}.region = i*ones(slength(x{i}),1);
  end
  x = concat_structs(x);
  return
end

if length(P.trailerseq)>P.trailerseq_len_required
  P.trailerseq = P.trailerseq(1:trailerseq_len_required);
end

% load and convert data
x = sffread(fname);
x = convert_sff(x);   % also suggests 'trim' for position at which to trim poor-quality bases
x = rename_field(x,'trim','qualtrim');
nx = slength(x);

% try to guess leader sequence
if nx>=100
  min_leader_len = 6;
  max_leader_len = 15;
  ltt = max_leader_len - min_leader_len + 1;
  tmp = randperm(nx);
  ntc = 100;
  check = tmp(1:ntc);
  z = cell(ntc,ltt);
  for i=1:ntc
    s = x.seq{check(i)};
    ls = length(s);
    for j=1:ltt
      z{i,j} = s(1:min(ls,min_leader_len-1+j));
    end
  end
  nu = nan(ltt,1);
  for j=1:ltt, nu(j) = length(unique(z(:,j))); end
  idx = find(nu==1,1,'last');
  imp = z{1,idx};
  fprintf('Imputed leader sequence: %s\n',imp);
  if ~isempty(P.leaderseq)
    if strcmp(P.leaderseq,imp)
      fprintf('(matches leader sequence specified in parameter struct)\n');
    else
      fprintf('using %s instead\n', P.leaderseq);
    end
  else
    fprintf('Will trim this leader from sequences\n');
    P.leaderseq = imp;
  end
end

% trim leaderseq
if isempty(P.leaderseq)
  x.leadertrim = zeros(nx,1);
else
  fprintf('Trimming leaders... ');
  x.leadertrim = nan(nx,1);
  ll = length(P.leaderseq);
  lim = ll+P.max_extra_trim;
  for i=1:nx, if ~mod(i,10000), fprintf('%d/%d ',i,nx); end
    sl = length(x.seq{i});
    fr = 1;
    to = min(sl,lim);
    p = strfind(x.seq{i}(fr:to),P.leaderseq);
    if length(p)>1
      if P.trim_extra_leaders, p = p(end); else p = p(1); end
    end
    if ~isempty(p), x.leadertrim(i) = p+ll-1; end
  end, fprintf('\n');
end

% trim trailerseq
if isempty(P.trailerseq)
  x.trailertrim = zeros(nx,1);
else
  fprintf('Trimming trailers... ');
  x.trailertrim = nan(nx,1);
  ll = length(P.trailerseq);
  lim = ll+P.max_extra_trim;
  for i=1:nx, if ~mod(i,10000), fprintf('%d/%d ',i,nx); end
    sl = length(x.seq{i});
    fr = max(1,x.qualtrim(i)-lim);
    to = min(sl,x.qualtrim(i)+10);
    p = strfind(x.seq{i}(fr:to),P.trailerseq);
    if length(p)>1
      if P.trim_extra_trailers, p = p(1); else p = p(end); end
    end
    if isempty(p)   % try smith-waterman
      [a b c] = swalign(x.seq{i}(fr:to),P.trailerseq);
      if a>=P.trailerseq_score_required, p=c(1); end
    end
    if ~isempty(p), x.trailertrim(i) = fr+p-1; end
  end, fprintf('\n');
end

% actually perform trimming
if ~P.soft_trim
  fprintf('Finalizing trimming: ');
  for i=1:nx, if ~mod(i,10000), fprintf('%d/%d ',i,nx); end
    if isnan(x.leadertrim(i)), lt=13; else lt=x.leadertrim(i); end
    x.seq{i} = x.seq{i}(lt+1:x.qualtrim(i));
    x.qual{i} = x.qual{i}(lt+1:x.qualtrim(i));
    x.len(i) = x.qualtrim(i)-lt;
  end, fprintf('\n');
  x = rmfield(x,{'qualtrim','leadertrim','trailertrim'});
  x = reorder_struct(x,x.len>=200);   % only keep sequences over a minimum length (after trimming)
end





