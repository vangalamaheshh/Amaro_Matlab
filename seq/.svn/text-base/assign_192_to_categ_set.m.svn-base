function c = assign_192_to_categ_set(K,C192)
% c = assign_192_to_categ_set(K,C192)
%
% Given a set of k categories as struct K with the following fields:
%   left   = subset of 'ACGT', representing 5' base
%   right  = subset of 'ACGT', representing 3' base
%   base   = subset of 'ACGT', representing mutated base
%            NOTE: if from contains only A&C, then will assume we are doing strand collapse
%            NOTE: "base" can also be called "from"
%   change = subset of 'tfs', representing Transition, Flip transversion, Skew transversion
%     OR
%   newbase= subset of 'ACGT' but not <base>, representing the base changed to
%            if length(base)>1, then newbase cannot be used.
%
% and the set of 192 elementary mutation types, as C192
%   left   = one of 'ACGT', representing 5' base
%   right  = one of 'ACGT', representing 3' base
%   base   = one of 'ACGT', representing mutated base
%   newbase= one of 'ACGT' but not <base>, representing the base changed to
% NOTE: if C192 is not provided, will assume the conventional default ordering.
%
% Maps C192 onto K.
%
% Returns vector c, with 192 rows, each of them 1-k, telling which of the k categories it belongs to
%
% NOTE: nulls/indels are not considered here.
% NOTE: gives an error if the categories in K are not surjective onto C192
%
% Mike Lawrence 2013-04-10

if nargin<2
  % assume default ordering of the 192 categories
  C192 = load_struct('/xchip/cga/reference/mutsig_params/192categs.txt');
end
% check C192
if isfield(C192,'from') && ~isfield(C192,'base'), C192=rename_field(C192,'from','base'); end
require_fields(C192,{'left','right','base','newbase'});
if isfield(C192,'change'), error('C192.change not supported'); end

% check K
k = slength(K);
if isfield(K,'from') && ~isfield(K,'base'), K=rename_field(K,'from','base'); end
require_fields(K,{'left','right','base'});
if ~isfield(K,'change') && ~isfield(K,'newbase'), error('need change or newbase'); end
if isfield(K,'type')
  pointidx = find(strcmpi(K.type,'point'));
else
  pointidx = 1:k;
end
for ki=1:k
  if length(K.base{ki})>1 && isfield(K,'newbase'), error('error in K'); end
end

% find out if K assumes strand collapse
u = unique(cat(2,K.base{:}));
if strcmp(u,'AC')
  % assume strand collapse
  strand_collapse = true;
elseif strcmp(u,'ACGT')
  % keep strands separate
  strand_collapse = false;
else
  error('problem with K.base');
end

base = 'ACGT';
complement(base) = 'TGCA';
whatchange('A','ACGT') = 'nstf';
whatchange('C','ACGT') = 'snft';
whatchange('G','ACGT') = 'tfns';
whatchange('T','ACGT') = 'ftsn';

% evaluate whether each category in K and each category in C192 are a match
m = false(192,k);
for ci=1:192
  base = C192.base{ci};
  newbase = C192.newbase{ci};
  left = C192.left{ci};
  right = C192.right{ci};
  if ismember(base,'GT') && strand_collapse==true
    base = complement(base);
    newbase = complement(newbase);
    left = complement(C192.right{ci});
    right = complement(C192.left{ci});
  end
  for ki=1:k
    % decide if the context belongs to this category
    if ~ismember(ki,pointidx), continue; end
    if ~ismember(base,K.base{ki}), continue; end
    if ~ismember(left,K.left{ki}), continue; end
    if ~ismember(right,K.right{ki}), continue; end
    % last match to verify: newbase
    match1 = nan;
    match2 = nan;
    if isfield(K,'newbase')
      match1 = 1*ismember(newbase,K.newbase{ki});
    end
    if isfield(K,'change')
      change = whatchange(base,newbase);
      match2 = 1*ismember(change,K.change{ki});
    end
    if (match1==0 && match2==1) || (match1==1 && match2==0)
      error('inconsistency between K.change and K.newbase');
    end
    if match1==1 | match2==1
      m(ci,ki)=1;
    end
  end
end

% make sure K->C192 mapping is surjective
sm = sum(m,2);
if any(sm>1)
  error('some categories in K map to more than one category in C192');
end
if any(sm<1)
  error('some categories in C192 are mapped to by no category in K');
end

% everything is good!
c = m*(1:k)';


%pr(C192.name,m,c,K.name(c))






