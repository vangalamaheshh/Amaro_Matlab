function F = find_double_nulls_in_fullM(F)

null_idx = find(strcmp('null',F.cat.name) | strcmp('indel+null',F.cat.name));
if isempty(null_idx), error('no null category found!'); end

double_null_idx = find(strcmp('double_null',F.cat.name) | strcmp('null2',F.cat.name));
if ~isempty(double_null_idx), error('double_null category already exists!'); end

% add double_null category to F.cat
F.cat = reorder_struct(F.cat,[1:slength(F.cat) slength(F.cat)]);
double_null_idx = slength(F.cat);
F.cat.autoname{double_null_idx} = 'double_null';
F.cat.name{double_null_idx} = 'double_null';

% remove per-category fields from F.ttype  (because now they will be confusing with 8 categories)
F.ttype = keep_fields_that_exist(F.ttype,{'name','npat','label','mutsig','short','npat_wgs'});

% find double-null mutations
mutflds = grep('^mut',fieldnames(F));
for i=1:length(mutflds)
  m = getfield(F,mutflds{i});
  m.is_double_null = false(slength(m),1);
  midx = find(m.categ==null_idx);
  [pu pi pj] = unique(m.patient(midx));
  [gu gi gj] = unique(m.gene(midx));
  [u ui uj] = unique([pj gj],'rows');
  h = histc(uj,1:length(u));
  ct = h(uj);
  idx = midx(ct>=2);
  m.is_double_null(idx) = true;
  m.categ(idx) = double_null_idx;
  F = setfield(F,mutflds{i},m);
end

