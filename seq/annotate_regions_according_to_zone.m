function R = annotate_regions_according_to_zone(R,P)
% R = annotate_regions_according_to_zone(R,P)
%
% annotates each region as to exonic/intronic/IGR/UTR

if ischar(P)
  tmp = P;
  P=[];
  P.build = tmp;
end

if ~exist('P','var'), P=[]; end

if ~isfield(P,'build')
  P.build = 'hg18';
  fprintf('Assuming hg18\n');
end

if strcmp(P.build,'hg18')
  dr = '/xchip/tcga_scratch/lawrence/db/zone';
  fldname = 'zone';
elseif strcmp(P.build,'hg19')
  dr = '/xchip/tcga_scratch/lawrence/db/hg19/zone';
  fldname = 'categ';
else
  error('unknown P.build %s',P.build);
end

P = impose_default_value(P,'ambiguity_resolution_method','max_overlap');

demand_fields(R,{'chr','start','end'});

if ~isnumeric(R.chr), R.chr = convert_chr(R.chr); end
if ~isnumeric(R.start), R.start = str2double(R.start); end
if ~isnumeric(R.end), R.end = str2double(R.end); end

nr = slength(R);

z = nan(nr,4);    % 0=IGR  1=intron  2=UTR  3=exon
for c=1:24, fprintf('%d/%d ',c,24);
  tmp = load([dr '/chr' num2str(c) '.mat']);
  zone = getfield(tmp,fldname);
  for i=1:nr
    if R.chr(i)==c
      st = min(length(zone),R.start(i));
      en = min(length(zone),R.end(i));
      z(i,:) = histc(zone(st:en),0:3);
end,end,end, fprintf('\n');

% resolve ambiguities
if strcmpi(P.ambiguity_resolution_method,'max_overlap')
  [tmp idx] = max(z,[],2);
elseif strcmpi(P.ambiguity_resolution_method,'zone_priority')
  for i=1:4, z(z(:,i)>0,i) = i; end
  [tmp idx] = max(z,[],2);
else
  error('unknown P.ambiguity_resolution_method %s',P.ambiguity_resolution_method);
end

zname = {'IGR';'intron';'UTR';'exon'};
R.zone = zname(idx);

% for regions in genes (exonic/intronic/UTR), find out what gene it is

D = load_refseq(P.build);
D.chrno = convert_chr(D.chr);
cidx = cell(24,1); for c=1:24, cidx{c} = find(D.chrno==c); end

R.gene = repmat({'-'},nr,1);
for i=1:nr, if ~mod(i,100000), fprintf('%d/%d ',i,nr); end
  if strcmp(R.zone{i},'IGR'), continue; end
  idx = cidx{R.chr(i)};
  if strcmp(R.zone{i},'UTR')
    idx = idx(D.tx_start(idx)-3000<=R.end(i));
    idx = idx(D.tx_end(idx)+3000>=R.start(i));
  else
    idx = idx(D.code_start(idx)<=R.end(i));
    idx = idx(D.code_end(idx)>=R.start(i));
    if strcmp(R.zone{i},'exon')
      keep = false(length(idx),1);
      for j=1:length(idx)
        e = find(D.exon_starts{idx(j)}<=R.end(i) & D.exon_ends{idx(j)}>=R.start(i));
        if ~isempty(e), keep(j)=true; break; end
      end
      idx = idx(keep);
    end
  end
  if isempty(idx), g = '?';
  elseif length(idx)==1, g = D.gene{idx};
  else g = concat(unique(D.gene(idx)),'/');
  end
  R.gene{i} = g;
end, fprintf('\n');


