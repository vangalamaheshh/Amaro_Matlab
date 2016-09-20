function M2 = unique_by_patient_and_site(M)
% M2 = unique_by_patient_and_site(M)
%
% Mike Lawrence 2010-01-27

fprintf('Uniquing by patient+site: ');

patient = M.patient;
chr = M.chr;
if isnumeric(chr), chr = chrlabel(chr); end
pos = M.start;
if isnumeric(pos), pos = num2cellstr(pos); end
site = stringsplice([chr,pos],1,':');
id = stringsplice([patient,site],1,'|');
[u ui uj] = unique(id);
M2 = reorder_struct(M,ui);

% concatenate certain fields
flds = {'dataset','sample_dir','center','Center'};
flds = intersect(flds,fieldnames(M2));
f = cell(length(flds),1);
for fi=1:length(flds), f{fi} = getfield(M,flds{fi}); end
f2 = cell(length(flds),1);;
for i=1:length(u)
  idx = find(uj==i);
  for fi=1:length(flds), f2{fi}{i,1} = concat(f{fi}(idx),';'); end
end
for fi=1:length(flds), M2 = setfield(M2,flds{fi},f2{fi}); end

fprintf('kept %d/%d records\n',slength(M2),slength(M));
