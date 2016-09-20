function X = make_genesets(in,genelist)

if ~exist('genelist','var')
  genelist = '/xchip/cga1/lawrence/capture/Refseq_exons_good_20101221_genelist.txt';
end

if ischar(in), in = {in}; end

g = load_struct(genelist);
demand_fields(g,'name');

X = [];
for i=length(in)
  q = load_lines(in{i});
  q = regexprep(q,'\s','');
  idx = listmap(q,g.name);
  fprintf('pathway %d:  %d/%d mapped\n',i,sum(~isnan(idx)),length(idx));
  if any(isnan(idx))
    fprintf('  MISSING:\n'); disp(q(isnan(idx)));
  end
  idx(isnan(idx)) = [];
  X.name{i,1} = ['pathway' num2str(i)];
  X.desc{i,1} = X.name{i};
  X.genes(i,1) = stringsplice(g.name(idx),2,char(9));
end

