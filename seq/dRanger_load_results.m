function X = dRanger_load_results(sample,P)
% dRanger_load_results(sample,P)

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'results_name',[]);

direc = ['/xchip/tcga_scratch/lawrence/' sample];
name2 = upper(regexprep(sample,'/','-'));
if isempty(P.results_name), results_name = [name2 '_dRanger_results'];
else results_name = P.results_name; end
fname = [direc '/' results_name '_bfiltered.txt'];
if ~exist(fname,'file'), error('Can''t find file %s',fname);end
X = load_struct(fname,[repmat('%f',1,14) repmat('%s',1,12) repmat('%f',1,25)]);

