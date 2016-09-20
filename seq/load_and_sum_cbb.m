function cov = load_and_sum_cbb(tumcbb,normcbb)
% loads each pair of tum+norm cbb's and applys the 14+8 cutoffs to get covered positions
% sums positions across genome, returns cell array of chromosomes

TUMCUTOFF = 14;
NORMCUTOFF = 8;

if length(tumcbb)~=length(normcbb), error('length(tumcbb)~=length(normcbb)'); end

for i=1:length(tumcbb)
  for chr=1:24
    tfile{i}{chr} = [tumcbb{i} '/chr' num2str(chr) '.txt'];
    nfile{i}{chr} = [normcbb{i} '/chr' num2str(chr) '.txt'];
  end
  demand_file([tfile{i};nfile{i}]);
end

cov = cell(24,1);
for i=1:length(tumcbb)
  fprintf('Individual %d/%d: ',i,length(tumcbb));
  for chr=1:24
    fprintf('chr%d ',chr);
    t = read_table(tfile{i}{chr},'%f',char(9),0);
    n = read_table(nfile{i}{chr},'%f',char(9),0);
    t = t.dat{1}; n = n.dat{1};
    if length(t)>length(n), n(length(t)) = 0; end
    if length(n)>length(t), t(length(n)) = 0; end
    c = double(t>=TUMCUTOFF & n>=NORMCUTOFF);
    if i==1
      cov{chr} = c;
    else
      if length(c)>length(cov{chr}), cov{chr}(length(c)) = 0; end
      if length(cov{chr})>length(c), c(length(cov{chr})) = 0; end
      cov{chr} = cov{chr} + c;
    end
  end
end













