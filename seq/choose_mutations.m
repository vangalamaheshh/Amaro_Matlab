function M = choose_mutations(M,P)

% default parameter values

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'use_all_mutations',false);
P=impose_default_value(P,'use_representative_of_repetitive',true);
P=impose_default_value(P,'use_germline_or_LOH',false);
P=impose_default_value(P,'use_repetitive',false);
P=impose_default_value(P,'use_indels',true);
P=impose_default_value(P,'use_validated_mutations_only',false);
P=impose_default_value(P,'use_validated_or_verified_mutations_only',false);
P=impose_default_value(P,'but_use_all_silent_mutations',false);

% choose mutations

use = find(~isnan(M.mut.patient) & ~isnan(M.mut.gene));

syn = grep('^Synonymous$', M.mut.type, 1);
mis = grep('^Missense$', M.mut.type, 1);
non = grep('^Nonsense$', M.mut.type, 1);
oth = grep('(In|Del|Frameshift|Splice|Targeted_Region)', M.mut.type, 1);

silent = syn;
nonsilent = union([mis;non],oth);

if ~P.use_all_mutations

  ok = find(strcmpi(M.mut.filtered,'OK'));
  rr = find(strcmpi(M.mut.filtered,'Representative_of_repetitive'));
  rep = find(strcmpi(M.mut.filtered,'Repetitive'));
  gol = find(strcmpi(M.mut.filtered,'Germline_or_LOH'));

  nonfiltered = ok;

  if P.use_representative_of_repetitive
     fprintf('Using representative_of_repetitive mutations\n');
     nonfiltered = union(nonfiltered, rr);
  end

  if P.use_repetitive
     fprintf('Using repetitive mutations\n');
     nonfiltered = union(nonfiltered,rep);
  end

  if P.use_germline_or_LOH
     fprintf('Using germline_or_LOH mutations\n');
     nonfiltered = union(nonfiltered,gol);
  end

  use = intersect(use, nonfiltered);

  if ~P.use_indels
    fprintf('Ignoring indels\n');
    indels = find(strcmpi(M.mut.class, 'Indel'));
    use = setdiff(use, indels);
  end

  use_tmp = use;

  if P.use_validated_mutations_only
    fprintf('Using validated mutations only\n');
    valid = grep('Validated', M.mut.evidence, 1);
    use = intersect(use, valid);
  end

  if P.use_validated_or_verified_mutations_only
    fprintf('Using validated/verified mutations only\n');
    vv = grep('Validated|Verified', M.mut.evidence, 1);
    use = intersect(use, vv);
  end

  if P.but_use_all_silent_mutations
    fprintf('But using all silent mutations\n');
    reinstate = intersect(use_tmp,silent);
    use = union(use,reinstate);
  end

else

  fprintf('Using all mutations\n');

end

M.use = use;
M.use_nonsilent = intersect(use,nonsilent);     % nonsilent
M.use_nonsense = intersect(use,non);            % nonsense
M.use_silent = intersect(use,silent);           % silent

end
