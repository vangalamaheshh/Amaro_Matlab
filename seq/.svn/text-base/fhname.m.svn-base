function f = fhname(samples)

if ~iscell(samples)
  samples = {samples};
  flag = true;
else
  flag = false;
end

for i=1:length(samples)
  name2 = upper(regexprep(samples{i},'/','-'));
  name2 = regexprep(name2,'-WGS','/wgs');
  name2 = regexprep(name2,'-CAPTURE','/capture');
  name2 = regexprep(name2,'MM-0028','MMRC0028');
  f{i,1} = name2;
end

if flag
  f = f{1};
end


