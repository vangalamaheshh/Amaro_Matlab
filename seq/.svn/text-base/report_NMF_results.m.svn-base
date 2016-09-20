function report_NMF_results(res)

if ~iscell(res)
  nf = size(res.h,1);
  tmp = res;
  clear res;
  res{nf} = tmp;
end

nfs = find(~cellfun('isempty',res));

for nf=nfs, fprintf('\n\n****************  %d FACTORS  ***********************\n',nf);
  for fi=1:nf, fprintf('\n%%**  FACTOR %d/%d  *****\n',fi,nf);;
    for k=2:5
      fprintf('\nmutation types (k=%d)\n\n',k);
      for i=1:slength(res{nf}.breakdown{fi}{k})
        fprintf('%20s  (relative rate %.2f)\n',res{nf}.breakdown{fi}{k}.name{i},res{nf}.breakdown{fi}{k}.relrate(i));
      end
    end
    fprintf('\ntumor types\n\n');
    for i=1:length(res{nf}.ttypefrac_text{fi});
      fprintf('%s\n',res{nf}.ttypefrac_text{fi}{i});
    end
  end
end

