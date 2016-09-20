function [G IGR] = analyze_intron_rates(X,N,G,Z,C,P)
% analyze_intron_rates(X,N,G,Z,C,P)
%
% input:
%
% X = mutation struct
% N = cell(24,1) of coverage by chromosome (number of covered samples at each basepair)
% G = list of genes (from get_refseq_introns_and_exons)
% Z = list of category names
% C = cell(24,1) of category by chromosome (category at each basepair)
% P = parameters struct
%
% output:
%
% G = with N and n added for introns and exons
% IGR = data for IGR
%
% Mike Lawrence 2010-06-20

try

if ~exist('P','var'), P=[]; end

require_fields(X,{'categnum','tum_allele1','tum_allele2','chr','start','end'});
require_fields(G,{'chr','footprint','n_exons','exon_starts','exon_ends','n_introns','intron_starts','intron_ends'});
if ~iscell(N) || length(N)~=24, error('N should be cell(24,1)'); end
if ~iscell(C) || length(C)~=24, error('C should be cell(24,1)'); end
require_fields(Z,{'name'});

maxsamps = double(max(cellfun(@max,N)));
ng = slength(G);
z = cell(ng,1);
G.exon_N = z; G.exon_n = z;
G.intron_N = z; G.intron_n = z;
G.tot_exon_N = z; G.tot_exon_n = z;
G.tot_intron_N = z; G.tot_intron_n = z;

% versions of analysis to do:
% (1) nothing excluded
% (2) exclude "bad" and "N": all territory: all mutations
% (3) exclude "bad" and "N": all territory: CpG transitions only
% (4) exclude "bad" and "N": all territory: mutations other than CpG transitions
% (5) exclude "bad" and "N": conserved territory only: all mutations
% (6) exclude "bad" and "N": conserved territory only: CpG transitions only
% (7) exclude "bad" and "N": conserved territory only: mutations other than CpG transitions
% (8) exclude "bad" and "N": nonconserved territory only: all mutations
% (9) exclude "bad" and "N": nonconserved territory only: CpG transitions only
%(10) exclude "bad" and "N": nonconserved territory only: mutations other than CpG transitions

% territory categories for each version
nver = 10;
zidx = cell(nver,1);
z = (1:slength(Z))'; zidx{1} = z;
z1 = z(grepv('bad|any N',Z.name(z),1)); zidx(2:4) = repmat({z1},3,1);
z1a = z1(grep('C in ._G|G in C_.',Z.name(z1),1)); zidx{3} = z1a;
z2 = z1(grepv('noncons',Z.name(z1),1)); zidx(5:7) = repmat({z2},3,1);
z2a = z2(grep('C in ._G|G in C_.',Z.name(z2),1)); zidx{6} = z2a;
z2 = z1(grep('noncons',Z.name(z1),1)); zidx(8:10) = repmat({z2},3,1);
z2a = z2(grep('C in ._G|G in C_.',Z.name(z2),1)); zidx{9} = z2a;

% mutation categories for each version
% use_all_mutations = 1;
% use_only_CpG_transitions = 2;
% use_all_except_CpG_transitions = 3;
cpgidx = [{[1 2]};repmat([{[1 2]};{2};{1}],3,1)];
X = make_numeric(X,{'categnum','chr','start','end'});
X.is_CpG_transition = false(slength(X),1);
X.is_CpG_transition(ismember(X.categnum,grep('C in ._G',Z.name,1))&(strcmp('T',X.tum_allele1)|strcmp('T',X.tum_allele2))) = true;
X.is_CpG_transition(ismember(X.categnum,grep('G in C_.',Z.name,1))&(strcmp('A',X.tum_allele1)|strcmp('A',X.tum_allele2))) = true;

% make chromosome lookup tables
cidx = cell(24,1); for i=1:24, cidx{i} = find(X.chr==i); end;

% get mutation counts and coverage counts for each intron and exon
fprintf('Analyzing genes... ');
for i=1:ng, if ~mod(i,100), fprintf('%d/%d ',i,ng); end
  if strcmp(G.name,'__IGR__'), continue; end
  chr = G.chr(i); chrlen = min(length(C{chr}),length(N{chr}));
  z = nan(nver,G.n_exons(i)); G.exon_N{i} = z; G.exon_n{i} = z;
  for j=1:G.n_exons(i)
    st = min(chrlen,G.exon_starts{i}(j));
    en = min(chrlen,G.exon_ends{i}(j));
    h = hist2d_fast(C{chr}(st:en),double(N{chr}(st:en)),1,slength(Z),0,maxsamps) * [0:maxsamps]';
    midx = cidx{chr}(X.start(cidx{chr})>=st & X.start(cidx{chr})<=en);
    mh = hist2d_fast(X.categnum(midx),double(X.is_CpG_transition(midx)),1,slength(Z),0,1);
    for v=1:nver
      G.exon_N{i}(v,j) = sum(h(zidx{v}));
      G.exon_n{i}(v,j) = sum(sum(mh(zidx{v},cpgidx{v})));
    end
  end
  G.tot_exon_N{i} = sum(G.exon_N{i},2);
  G.tot_exon_n{i} = sum(G.exon_n{i},2);
  z = nan(nver,G.n_introns(i)); G.intron_N{i} = z; G.intron_n{i} = z;
  for j=1:G.n_introns(i)
    st = min(chrlen,G.intron_starts{i}(j));
    en = min(chrlen,G.intron_ends{i}(j));
    h = hist2d_fast(C{chr}(st:en),double(N{chr}(st:en)),1,slength(Z),0,maxsamps) * [0:maxsamps]';
    midx = cidx{chr}(X.start(cidx{chr})>=st & X.start(cidx{chr})<=en);
    mh = hist2d_fast(X.categnum(midx),double(X.is_CpG_transition(midx)),1,slength(Z),0,1);
    for v=1:nver
      G.intron_N{i}(v,j) = sum(h(zidx{v}));
      G.intron_n{i}(v,j) = sum(sum(mh(zidx{v},cpgidx{v})));
    end
  end
  G.tot_intron_N{i} = sum(G.intron_N{i},2);
  G.tot_intron_n{i} = sum(G.intron_n{i},2);
end, fprintf('\n');

% get total n and N for IGR
fprintf('Computing IGR rate... ');
IGR = [];
IGR.tot_n = nan(nver,1);
mh = hist2d_fast(X.categnum,double(X.is_CpG_transition),1,slength(Z),0,1);
for v=1:nver
  zidx_igr{v,1} = zidx{v}(grep('IGR',Z.name(zidx{v}),1));
  IGR.tot_n(v) = sum(sum(mh(zidx_igr{v},cpgidx{v})));
end
IGR.tot_N = zeros(nver,1);
for chr=1:24, fprintf('chr%d ',chr);
  minlen = min(length(C{chr}),length(N{chr}));
  h = hist2d_fast(C{chr}(1:minlen),double(N{chr}(1:minlen)),1,slength(Z),0,maxsamps) * [0:maxsamps]';
  for v=1:nver, IGR.tot_N(v) = IGR.tot_N(v) + sum(h(zidx_igr{v})); end
end, fprintf('\n');

fprintf('DONE\n');

catch me; excuse(me); end
