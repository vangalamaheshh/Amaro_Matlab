function M = collapse_multiM(Ms)

check_multiM_agreement(Ms);
nm = length(Ms);

M = [];
M.name = 'combined';

% report
M.report = [];
for i=1:nm
  M.report = [M.report Ms{i}.report];
end

% cov
M.cov = [];
tmp = cell(nm,1);
for i=1:nm
  tmp{i} = Ms{i}.cov.sample;
  tmp{i}.isetname = repmat({Ms{i}.name},slength(tmp{i}),1);
end
tmp = concat_structs_keep_all_fields(tmp);
if length(unique(tmp))<length(tmp), error('Non-unique patient names in cov!'); end
M.cov.sample = tmp;
M.cov.ns = slength(M.cov.sample);
M.cov.targ = Ms{1}.cov.targ;
M.cov.nt = Ms{1}.cov.nt;
tmp = cell(nm,1); for i=1:nm, tmp{i} = Ms{i}.cov.fcov; end 
M.cov.fcov = cat(2,tmp{:});
if isfield(Ms{1}.cov,'gene')
  M.cov.gene = Ms{1}.cov.gene;
  M.cov.ng = Ms{1}.cov.ng;
  tmp = cell(nm,1); for i=1:nm, tmp{i} = Ms{i}.cov.gene_totcov; end
end
M.cov.gene_totcov = cat(2,tmp{:});

% gene
M.gene = Ms{1}.gene;
M.ng = Ms{1}.ng;

% patient
tmp = cell(nm,1);
for i=1:nm
  tmp{i} = Ms{i}.patient;
  tmp{i}.isetname = repmat({Ms{i}.name},slength(tmp{i}),1);
end
tmp = concat_structs_keep_all_fields(tmp);
if length(unique(tmp))<length(tmp), error('Non-unique patient names in cov!'); end
M.patient = tmp;
M.np = slength(M.patient);
M.patient.cidx = listmap(M.patient.name,M.cov.sample.name);
if isfield(M.patient,'pidx'), M.patient = rmfield(M.patient,'pidx'); end

% mut, use, use_nonsilent
tmp = cell(nm,1);
usesil = cell(nm,1);
usenon = cell(nm,1);
offset = 0;
for i=1:nm
  tmp{i} = Ms{i}.mut;
  tmp{i}.isetname = repmat({Ms{i}.name},slength(tmp{i}),1);
  usesil{i} = Ms{i}.use_silent + offset;
  usenon{i} = Ms{i}.use_nonsilent + offset;
  offset = offset + slength(Ms{i}.mut);
end
tmp = concat_structs_keep_all_fields(tmp);
M.mut = tmp;
M.mut.patient = listmap(M.mut.patient_name,M.patient.name);
M.use_silent = cat(2,usesil{:});
M.use_nonsilent = cat(2,usenon{:});
M.use = union(M.use_silent,M.use_nonsilent);

% N_cov, n_silent, n_nonsilent
M.TOT = 1;
M.mutclass = {'Total'};
M.NUM_INDEL_CLASSES = 0;
M.N_terr = Ms{1}.N_terr(:,Ms{1}.TOT);
N = cell(nm,1);
nsil = cell(nm,1);
nnon = cell(nm,1);
for i=1:nm
  N{i} = Ms{i}.N_cov(:,Ms{i}.TOT,:);
  nsil{i} = Ms{i}.n_silent(:,Ms{i}.TOT,:);
  nnon{i} = Ms{i}.n_nonsilent(:,Ms{i}.TOT,:);
end
M.N_cov = cat(3,N{:});
M.n_silent = cat(3,nsil{:});
M.n_nonsilent = cat(3,nnon{:});

