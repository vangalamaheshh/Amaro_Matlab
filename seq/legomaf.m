function [z c] = legomaf(X,P)
% legomaf(X,P)

if exist('P','var') && ischar(P)
  fprintf('Using "%s" as the coverage model.\n',P);
  tmp=P;
  P=[];
  P.coverage = tmp;
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'display','rates');
P = impose_default_value(P,'coverage','exome');

if ~isfield(X,'newbase'), X.newbase = find_newbase(X); end

require_fields(X,{'context65','newbase'});
X = make_numeric(X,'context65');
X.newbase_idx = listmap(X.newbase,{'A','C','G','T'});

if isfield(X,'classification'), X = reorder_struct(X,strcmp('SNP',X.classification)); end
if isfield(X,'Variant_Type'), X = reorder_struct(X,strcmp('SNP',X.Variant_Type)); end
if isfield(X,'is_indel'), X = reorder_struct(X,~X.is_indel); end

n = hist2d_fast(X.context65,X.newbase_idx,1,64,1,4);

if grepmi('exome',{P.coverage}), N = average_exome_coverage();
elseif grepmi('genome|wgs',{P.coverage}), N = average_genome_coverage();
elseif grepmi('^(flat|counts|n)$',{P.coverage}), N = 1e6*ones(64,1);
else error('Unknown P.coverage %s',{P.coverage}); end

if isfield(X,'patient'), npat = length(unique(X.patient));
elseif isfield(X,'pat_idx'), npat = length(unique(X.pat_idx));
elseif isfield(X,'patient_idx'), npat = length(unique(X.patient_idx));
elseif isfield(X,'Tumor_Sample_Barcode'), npat = length(unique(X.Tumor_Sample_Barcode));
else npat=1; end

N=N*npat;

Nn = [N n];
Nn = collapse_Nn_64_by_strand(Nn);

% convert to 96-row format:
% 17-32 >A
% 1-16 >C
% 1-16 >G
% 17-32 >G
% 1-16 >T
% 17-32 >T
N = [Nn(17:32,1);Nn(1:16,1);Nn(1:16,1);Nn(17:32,1);Nn(1:16,1);Nn(17:32,1)];
n = [Nn(17:32,2);Nn(1:16,3);Nn(1:16,4);Nn(17:32,4);Nn(1:16,5);Nn(17:32,5)];

if strcmp(P.display,'counts'), data = n;
elseif strcmp(P.display,'rates'), data = n./N;
else error('Unknown P.display %s',P.display); end

% Lego plot
z = factor2lego(data');
lego(z,P);

if nargout==0
  clear z c
end
