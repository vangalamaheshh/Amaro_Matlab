function [x,report] = preprocess_mutations(x,categdir,P)

if nargin<2, error('requires x and categdir'); end

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'build','*required*');
P=impose_default_value(P,'build_dir','*required*');
P=impose_default_value(P,'consolidate_adjacent_muts',true);
P=impose_default_value(P,'consolidate_adjacent_muts_threshold',0);
P=impose_default_value(P,'remove_noncoding_mutations',true);
if ~isfield(P,'enforce_target_list') && isfield(P,'enforce_targlist')
  P.enforce_target_list = P.enforce_targlist;
end
P=impose_default_value(P,'enforce_target_list',false);
if P.enforce_target_list
  if ~isfield(P,'target_list') && isfield(P,'targlist')
    P.target_list = P.targlist;
  end
  if ~isfield(P,'target_list') && isfield(P,'targetlist')
    P.target_list = P.targetlist;
  end
  P=impose_default_value(P,'target_list','*required*');
end
P = impose_default_value(P,'mutation_blacklist',[]);
P = impose_default_value(P,'mutation_blacklist_match_fields','chr,start,newbase');
P = impose_default_value(P,'mutation_whitelist',[]);
P = impose_default_value(P,'mutation_whitelist_match_fields','patient,chr,start,newbase');

report = sprintf('Total number of mutations in input MAFs:\t%d\n',slength(x));

% copy some fields, add newbase, and convert some fields to numeric
x = add_and_convert_simple_fieldnames(x,P);

% remove invalidated mutations
if isfield(x,'Validation_Status')
  fprintf('Validation_Status:\n');
  count(x.Validation_Status);
  idx = grepi('Reference|LOH|Germline|Wild.?type|Artifact|Invalid|Bad|Remove|Do.?N.t.?Use',x.Validation_Status,1);
  if ~isempty(idx)
    fprintf('Removing the following invalidated mutations.');
    count(x.Validation_Status(idx));
    x = reorder_struct_exclude(x,idx);
    report = [report sprintf('After removing %d invalidated mutations:\t%d\n',length(idx),slength(x))];
  end
else
  fprintf('No Validation_Status column: unable to exclude any mutations on basis of validation data.\n');
end

% remove mutations outside the canonical chromosomes
ct = get_chrct(P.build);
idx = find(isnan(x.chr));
if ~isempty(idx)
  fprintf('Omitting the following %d mutations on nonstandard chromosomes:\n',length(idx));
  if isfield(x,'Chromosome'), count(x.Chromosome(idx)); else look(x,idx); end
  x = reorder_struct_exclude(x,idx);
  report = [report sprintf('After removing %d mutations outside chr1-%d:\t%d\n',length(idx),ct,slength(x))];
end

% remove mutations with non-numeric start/end
idx = find(isnan(x.start) | isnan(x.end));
if ~isempty(idx)
  fprintf('Omitting the following %d mutations with malformed start/end coordinates:\n',length(idx));
  look(x,idx);
  x = reorder_struct_exclude(x,idx);
  report = [report sprintf('After removing %d mutations with malformed start/end coordinates:\t%d\n',length(idx),slength(x))];
end

% apply mutation_blacklist (if available)
if ~isempty(P.mutation_blacklist)
  if ~exist(P.mutation_blacklist,'file')
    fprintf('Mutation blacklist file not found:  %s\n',P.mutation_blacklist);
  else
    try
      fprintf('Applying mutation blacklist:  %s\n',P.mutation_blacklist);
      B = load_struct(P.mutation_blacklist);
      if ~isfield(B,'patient'), B.patient = repmat({'null'},slength(B),1); end
      B = add_and_convert_simple_fieldnames(B);
      flds = split(P.mutation_blacklist_match_fields,',');
      if ~all(isfield(x,flds)) || ~all(isfield(B,flds))
        fprintf('Missing some fields specified in P.mutation_blacklist_match_fields = %s\n',flds);
      else
        bi = multimap(x,B,flds);
        idx = find(~isnan(bi));
        if ~isempty(idx)
          fprintf('Omitting the following %d blacklisted mutations:\n', length(idx));
          disp([x.patient(idx) x.gene(idx) num2cell([x.chr(idx) x.start(idx)]) x.newbase(idx)])
          x = reorder_struct_exclude(x,idx);
          report = [report sprintf('After removing %d blacklisted mutations:\t%d\n',length(idx),slength(x))];
        end
      end
    catch me
      fprintf('Error while attempting to apply blacklist:\n');
      disp(me.message);
      disp(me);
    end   
  end
end

% apply mutation_whitelist (if available)
if ~isempty(P.mutation_whitelist)
  if ~exist(P.mutation_whitelist,'file')
    fprintf('Mutation whitelist file not found:  %s\n',P.mutation_whitelist);
  else
    try
      fprintf('Applying mutation whitelist:  %s\n',P.mutation_whitelist);
      W = load_struct(P.mutation_whitelist);
      W = add_and_convert_simple_fieldnames(W);
      % first see if any mutations have to be overwritten
      flds = split(P.mutation_whitelist_match_fields,',');
      if ~all(isfield(x,flds)) || ~all(isfield(W,flds))
        fprintf('Missing some fields specified in P.mutation_whitelist_match_fields = %s\n',flds);
      else
        wi = multimap(x,W,flds);
        idx = find(~isnan(wi));
        if ~isempty(idx)
          fprintf('Overwriting the following %d whitelisted mutations:\n', length(idx));
          disp([x.patient(idx) x.gene(idx) num2cell([x.chr(idx) x.start(idx)]) x.newbase(idx)])
          x = reorder_struct_exclude(x,idx);
          report = [report sprintf('After removing %d mutations overwritten by whitelist:\t%d\n',length(idx),slength(x))];
        end
      end
      % then add the mutations from the whitelist
      fprintf('Adding the following %d whitelisted mutations:\n', slength(W));
      disp([W.patient W.gene num2cell([W.chr W.start]) W.newbase])
      x = concat_structs_keep_all_fields({x,W});
      report = [report sprintf('After adding %d mutations from whitelist:\t%d\n',slength(W),slength(x))];
    catch me
      fprintf('Error while attempting to apply whitelist:\n');
      disp(me.message);
      disp(me);
    end
  end
end

% ANNOTATE FLANK, CODING, SILENT, INDEL
x = add_helper_is_fields(x);

% REMOVE NONCODING MUTATIONS (if specified)
if P.remove_noncoding_mutations
  idx = find(~x.is_coding);
  if ~isempty(idx)
    fprintf('Removing the following non-coding mutations (based on "type" field):');
    count(x.type(idx));
    x = reorder_struct_exclude(x,idx);
    report = [report sprintf('After removing %d noncoding mutations:\t%d\n',length(idx),slength(x))];
   end
else
  fprintf('NOT removing non-coding mutations.\n');
end

% ENFORCE TARGET LIST
if P.enforce_target_list
  fprintf('Enforcing target list:\n',P.target_list);
  T = load_target_file(P.target_list,P);
  map = map_mutations_to_targets(x,T,P);
  idx = find(isnan(map));
  if ~isempty(idx)
    fprintf('Removing %d/%d mutations that fall outside target intervals.\n',length(idx),slength(x));
    x = reorder_struct_exclude(x,idx);
    map(idx)=[];
    report = [report sprintf('After removing %d mutations outside target intervals:\t%d\n',length(idx),slength(x))];
  end
  old_gene = x.gene;
  x.gene = T.gene(map);
  idx = find(~strcmpi(old_gene,x.gene));
  if ~isempty(idx)
    fprintf('Reassigning the following %d gene identities:', length(idx));
    tmp = stringsplice([old_gene(idx) x.gene(idx)],1,' -> ');
    count(tmp,1);
  end
else
  fprintf('Not enforcing any target list.\n');
end

% collapse adjacent mutations
if P.consolidate_adjacent_muts
  tmp = slength(x);
  x = collapse_adjacent_mutations(x,P);
  if slength(x)<tmp
    report = [report sprintf('After collapsing adjacent/redundant mutations:\t%d\n',slength(x))];
  end
end

% add context
x.context_orig = get_context(x.chr,x.start,categdir,P);

% collapse to context65 (if necessary)
categs_txt = [categdir '/categs.txt'];
try
  map65 = map_categories_to_65(categs_txt);
  x.context65 = nansub(map65,x.context_orig);
catch me
  disp(me); disp(me.message);
  error('ERROR when trying to collapse %s', categs_txt);
end

