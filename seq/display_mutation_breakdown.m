function display_mutation_breakdown(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'omit_silent',false);
P=impose_default_value(P,'mutation_breakdown_output',[]);

fprintf('Mutation breakdown by type\n');

if isfield(M,'use')
  q = M.mut.type(M.use);
else
  q = M.mut.type;
end

if P.omit_silent
  q=q(setdiff(1:length(q),grep('Syn',q,1)));
end


%q = regexprep(q,'(Frameshift|Inframe)_(Ins|Del)','Indel');
q = regexprep(q,'^Missense$','Missense_Mutation');
q = regexprep(q,'^Nonsense$','Nonsense_Mutation');

count(q);

if ~isempty(P.mutation_breakdown_output)
  [u ui uj] = unique(q);
  h = histc(uj,1:length(u));
  X = [];
  X.type = u;
  X.count = h;
  X.type{end+1,1} = 'Total';
  X.count(end+1,1) = sum(X.count);
  save_struct(X,P.mutation_breakdown_output);
end
