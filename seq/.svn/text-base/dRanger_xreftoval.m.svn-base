function X = dRanger_xreftoval(X,P)
% Mike Lawrence 2009-04-07
%
% rewritten 2010-03-31

if ~exist('P','var'), P=[]; end

if ischar(P), error('Second parameter should be P struct'); end


alias.old = {'OV-04-1371';'OV-13-0725';'OV-13-0751';'OV-24-0982';'OV-25-1319';...
   'GBM-06-0145';'GBM-06-0152';'GBM-06-0185';'GBM-06-0188';'GBM-06-0648';...
   'JN_08';'JN_97';'JN_TT';'MM-0421-FIX';'MM-0422-FIX';'MM-0425-FIX'};
alias.new = {'OV-1371';'OV-0725';'OV-0751';'OV-0982';'OV-1319';...
   'GBM-0145';'GBM-0152';'GBM-0185';'GBM-0188';'GBM-0648';...
   'CLL-JN';'CLL-JN';'CLL-JN';'MM-0421';'MM-0422';'MM-0425'};

P = impose_default_value(P,'method',2);
P = impose_default_value(P,'match_margin',400);
P = impose_default_value(P,'alias',alias);
P = impose_default_value(P,'valfile','/xchip/cga1/lawrence/dRanger/20091016_val-pr/all_validation_summary_20100424.txt');

if P.method==1
  % April 2009 method
  valfile='/xchip/tcga/gbm/analysis/lawrence/dRanger/dRanger_validation.txt';
  V = load_struct(valfile,'%s%f%f%f%f%f%f%f%f%s%s%s%s%s%s%s%s%s%s%s%s');

  if ~isempty(P.alias), error('P.alias not supported with P.method=1'); end
  
  nx = slength(X);
  z = cell(nx,1);
  X.name=z; X.well=z; X.A=z; X.B=z; X.R=z; X.val=z;
  
  for i=1:nx
    idx = find(   ( V.chr1==X.chr1(i) & V.min1<X.max1(i) & V.max1>X.min1(i) & V.str1==X.str1(i)    &...
                    V.chr2==X.chr2(i) & V.min2<X.max2(i) & V.max2>X.min2(i) & V.str2==X.str2(i) )   ...
                  | ( V.chr2==X.chr1(i) & V.min2<X.max1(i) & V.max2>X.min1(i) & V.str2==X.str1(i)    &...
                      V.chr1==X.chr2(i) & V.min1<X.max2(i) & V.max1>X.min2(i) & V.str1==X.str2(i) )       );
    X.name{i} = concat(V.id(idx),'/');
    X.well{i} = concat(V.well(idx),'/');
    X.A{i} = concat(V.A(idx),'/');
    X.B{i} = concat(V.B(idx),'/');
    X.R{i} = concat(V.R(idx),'/');
    X.val{i} = concat(V.val(idx),'/');
  end

elseif P.method==2
  % March 2010 method
  flds = {'individual','chr1','str1','pos1','chr2','str2','pos2'};
  require_fields(X,flds);
  X = make_numeric(X,setdiff(flds,{'individual'}));

  V = load_struct(P.valfile);
  require_fields(V,flds);
  V = make_numeric(V,setdiff(flds,{'individual'}));

  if ~isempty(P.alias)
    Vindividual = apply_aliases(V.individual,P.alias);
    Xindividual = apply_aliases(X.individual,P.alias);
  else
    Vindividual = V.individual;
    Xindividual = X.individual;
  end

  vidx = cell(slength(X),1);
  [u ui uj] = unique(Vindividual);
  mi = {};
  for i=1:slength(X)
    % match individual
    ii = find(strcmp(Xindividual{i},u));
    if isempty(ii)
      mi = [mi;X.individual{i}]; 
      continue;
    end
    idx = find(uj==ii);
    % match chr+str
    idx = idx(V.chr1(idx)==X.chr1(i));
    idx = idx(V.chr2(idx)==X.chr2(i));
    idx = idx(V.str1(idx)==X.str1(i));
    idx = idx(V.str2(idx)==X.str2(i));
    % match pos
    idx = idx(abs(V.pos1(idx)-X.pos1(i))<=P.match_margin);
    idx = idx(abs(V.pos2(idx)-X.pos2(i))<=P.match_margin);
    vidx{i} = idx;
  end

  if ~isempty(mi)
    fprintf('The database lacks validation info for the following individuals:\n');
    disp(unique(mi));
  end

  addflds = fieldnames(V);
  addflds(ismember(addflds,flds)) = [];
  vv = cell(slength(V),length(addflds));
  for j=1:length(addflds)
    vv(:,j) = getfield(V,addflds{j});
  end

  xx = repmat({'-'},slength(X),length(addflds));
  for i=1:slength(X)
    idx = vidx{i};
    if length(idx)==1
      xx(i,:) = vv(idx,:);
    elseif length(idx)>1
      for j=1:length(addflds)
        xx{i,j} = concat(vv(idx,j),';');
      end
    end
  end

  for j=1:length(addflds)
    X = setfield(X,addflds{j},xx(:,j));
  end

else
  error('Unknown P.method');
end
