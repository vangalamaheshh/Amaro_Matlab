function demand_multiP_files(P)

if ~iscell(P), error('P should be cell'); end

list = {};
for i=1:length(P)
  demand_field(P{i},'covfile'), list = [list;P{i}.covfile];
  demand_field(P{i},'mutfile'), list = [list;P{i}.mutfile];
  if isfield(P{i},'catfile'), list = [list;P{i}.catfile]; end
  if isfield(P{i},'patlist'), list = [list;P{i}.patlist]; end
  if isfield(P{i},'genelist'), list = [list;P{i}.genelist]; end
  if isfield(P{i},'targlist'), list = [list;P{i}.targlist]; end
  if isfield(P{i},'summed_cov_track'), list = [list;P{i}.summed_cov_track]; end
end

demand_file(list);
