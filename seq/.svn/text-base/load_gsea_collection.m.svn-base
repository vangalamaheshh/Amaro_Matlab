function C = load_gsea_collection(filename, remove_slashes)
%
% load_gsea_collection(filename, remove_slashes)
%
% Reads a Gene Set Enrichment Analysis collection from the specified file.
%
% File format is as follows:
%   No header line.
%   Each record line is:
%     Geneset_name Geneset_descrption Gene_1 Gene_2 Gene_3 ... Gene_n
%       (no predetermined number of genes on each line)
%
% Returns a structure C with the following fields:
%   name = Geneset_name
%   desc = Geneset_description
%   genes = cell array of gene names
%
% If parameter remove_slashes is used, then gene names of the following form:
%    abc /// xyz
%    abc_///_xyz
%    abc /// xyz /// ijk
% are converted to individual gene names (abc, xyz, ijk), and slashes are discarded.
%
% Mike Lawrence 2008-06-17
%

% read whole file in one gulp

in = fopen(filename);
F = fread(in,'char=>char')';
fclose(in);

% process file into genesets

S = split(F, char(10));   % split into lines
S = S(~cellfun('isempty',S));    % remove blank lines
nS = length(S);

C = [];
C.name = cell(nS,1);
C.desc = cell(nS,1);
C.genes = cell(nS,1);

% process each geneset

for s=1:nS
  if exist('remove_slashes','var')
    S{s} = regexprep(S{s},' /// |_///_',char(9));
  end
  fld = split(S{s}, char(9));  % split at tabs
  C.name{s} = fld{1};
  C.desc{s} = fld{2};
  C.genes{s} = fld(3:end); 
end
