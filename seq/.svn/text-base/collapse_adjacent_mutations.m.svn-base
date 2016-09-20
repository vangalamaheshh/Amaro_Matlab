function M = collapse_adjacent_mutations(M,P)
% M = collapse_adjacent_mutations(M,P)
% collapse adjacent mutations to single mutations
%
% Mike Lawrence 2010-01-27

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'use_new_simplified_procedure',1);
P=impose_default_value(P,'consolidate_adjacent_muts_threshold',1);

if P.use_new_simplified_procedure
%  fprintf('collapse_adjacent_mutations: using new simplified procedure.\n');
  fprintf('Will only collapse redundant mutations; will leave adjacent mutations as is.\n');

  [u ui uj] = unique_combos(M.patient,M.chr,M.start);
  if length(u)==slength(M)
    fprintf('All mutations are already unique: no collapsing required.\n');
    return
  end
  
  fprintf('Only %d/%d unique mutations: will collapse redundant mutations.\n',length(u),slength(M));
  keep = true(slength(M),1);
  for i=1:length(u)
    idx = find(uj==i);
    if length(idx)==1, continue; end
    z = unique(M.classification(idx));
    if length(z)==1
      keepnum = 1;
    else
      z = find(grepm('DEL|INS',M.classification(idx)),1);  % preferentially retain indels
      if isempty(z)
        keepnum=1;
      else
        keepnum=z;
      end
    end
    keep(idx) = false;
    keep(idx(keepnum)) = true;
  end

  M = reorder_struct(M,keep);
  fprintf('Done.\n');
  return
end

%%%%%%%%%%%%%%% OLDER PROCEDURE:  CAUSED PROBLEMS


if P.consolidate_adjacent_muts_threshold>=2
  fprintf('WARNING:  P.consolidate_adjacent_muts_threshold>=2; access to reference sequence is required.\n');
  if ~isfield(P,'build')
    fprintf('   Because reference sequence directory is not specified in P.build, default hg18 location is being used.\n');
  else
    if P.build(1)=='/'
      fprintf('   Using %s as reference sequence directory\n',P.build);
    else
      fprintf('   Using default reference sequence directory based on P.build=%s\n',P.build);
    end
  end
end
P=impose_default_value(P,'build','hg18');

require_fields(M,'classification');

if ~isfield(M,'ref_allele'), M.ref_allele = M.Reference_Allele; end
if ~isfield(M,'newbase'), M.newbase = find_newbase(M.newbase); end

fprintf('Collapsing adjacent mutations:\n');

numstart = slength(M);

done_collapsing = false;
while(~done_collapsing)

  n = nan(slength(M),1);
  idx = cell(slength(M),1);
  events = {}; genes = {}; patients = {};
  
  pat = M.patient;
  chr = M.chr;
  if isnumeric(chr), chr = chrlabel(chr); end
  M = make_numeric(M,{'start','end'});
  pos = M.start;
  
  pc = stringsplice([pat,chr],1,'|');
  [u ui uj] = unique(pc);
  for i=1:length(u)
    j = find(uj==i);
    x = dist(pos(j),[]);
    x = (x<=P.consolidate_adjacent_muts_threshold);
    for jj=1:length(j)
      c = j(x(:,jj));
      idx{j(jj)} = c;
      n(j(jj)) = length(c);
    end
  end
  
  keep = true(slength(M),1);
  M2 = M;

  done_collapsing = true;  
  for i=1:slength(M2), j=idx{i};
    if ~keep(i), continue; end
    if n(i)>1
      done_collapsing = false;
      M2.start(i) = min(M.start(j));
      M2.end(i) = max(M.end(j));
      len = M2.end(i) - M2.start(i) + 1;
      if any(strcmpi(M2.classification(j),'Ins')|strcmpi(M2.classification(j),'Del'))
        jdx = grepi('Ins|Del',M2.classification(j),1);
        newclass = M2.classification{jdx(1)};
        newtype = M2.type{jdx(1)};
        M2.classification{i} = newclass;
        M2.type{i} = newtype;
        event = sprintf('%d SNPs/Indels to one %s',length(j),newclass);
      else
        if len==1
          M2.classification{i} = 'SNP';
          event = sprintf('%d redundant SNPs to one SNP',length(j));
        elseif len==2
          M2.classification{i} = 'DNP';
          M2.type{i} = 'Missense';  % (approximation: need to actually compute this)
          event = sprintf('%d SNPs to one DNP (guessing "Missense")',length(j));
        else
          M2.classification{i} = 'Complex_substitution';
          M2.type{i} = 'Missense';  % (approximation: need to actually compute this)
          event = sprintf('%d SNPs to one Complex_substitution (guessing "Missense")',length(j));
        end
        if isfield(M2,'dataset')
          M2.dataset{i} = concat(unique(M2.dataset(j)),'+');
        end
      end
      events = [events;event];
      genes = [genes;M2.gene{i}];
      patients = [patients;M2.patient{i}];
      if P.consolidate_adjacent_muts_threshold < 2
        % we're only collapsing mutations at identical or adjacent sites: don't need access to reference sequence
        ref = repmat('N',1,M2.end(i)-M2.start(i)+1);
        for k=1:length(j), q = j(k);
          try
            ref(M.start(q)-M2.start(i)+1:M.end(q)-M2.start(i)+1) = M.ref_allele{q};
          catch me
            fprintf('Error processing mutation (ref sequence length not equal to (end-start+1):\n'); look(M,q);
          end
        end
      else
        % we're allowing collapse of two SNPs separated by a reference base: need access to reference sequence
        ref = upper(genome_region(M.chr(i),M2.start(i),M2.end(i),P.build));
      end
%      if isfield(M,'tum_allele1'), tum1_is_ref = strcmp(M.tum_allele1{i},M.ref_allele{i}); end
%      if isfield(M,'tum_allele2'), tum2_is_ref = strcmp(M.tum_allele2{i},M.ref_allele{i});
      M2.ref_allele{i} = ref;
      tum = ref;
      for k=1:length(j); q = j(k);
        newseq = M.newbase{q};
%        if ~strcmp(M.tum_allele1{q},M.ref_allele{q}), newseq = M.tum_allele1{q};
%        else newseq = M.tum_allele2{q};
%        end
        try
          tum = [tum(1:(M.start(q)-M2.start(i)))...
                 newseq...
                 tum(M.end(q)-M2.start(i)+2:end)];
          %tum([M.start(q):M.end(q)]-M2.start(i)+1) = newseq;
        catch me
          fprintf('Problem in collapse_adjacent_mutations\n');
        end
      end
%      if ~tum1_is_ref, M2.tum_allele1{i} = tum; else M2.tum_allele1{i} = ref; end
%      if ~tum2_is_ref, M2.tum_allele2{i} = tum; else M2.tum_allele2{i} = ref; end
      if isfield(M,'tum_allele1') && isfield(M,'tum_allele2')
        M2.tum_allele1{i} = tum;
        M2.tum_allele1{i} = ref;
      end
      if isfield(M,'newbase')
        M2.newbase{i} = tum;
      end
      if isfield(M,'genomechange')       
        if iscell(M2.chr)
          cc = chrlabel(M2.chr{i});
        else
          cc = chrlabel(M2.chr(i)); cc=cc{1};
        end
        M2.genomechange{i} = ['g.' cc ':' num2str(M2.start(i)) ref '>' tum];
        if iscell(M2.genomechange{i}), error('String/Cell concatenation problem'); end
        M2.cDNAchange{i} = '?';
        M2.codonchange{i} = '?';
        M2.proteinchange{i} = '?';
        end
        keep(j(2:end)) = false;
      end
  end
  
  M = reorder_struct(M2,keep);
  if ~isempty(events)
    count(events);
    fprintf('\nGenes affected:'); count(genes,1);
    fprintf('\nPatients affected:'); count(patients,1);
  end
  
end   % until done collapsing


numend = slength(M);

fprintf('Collapsed %d muts to %d muts\n',numstart,numend);

% copy the changes back into the original fields (if applicable)  
if isfield(x,'Variant_Classification') && isfield(x,'type'), x.Variant_Classification = x.type; end
if isfield(x,'Variant_Type') && isfield(x,'classification'), x.Variant_Type = x.classification; end
if isfield(x,'Start_position') && isfield(x,'start'), x.Start_position = x.start; end
if isfield(x,'End_position') && isfield(x,'end'), x.End_position = x.end; end
if isfield(x,'Reference_Allele') && isfield(x,'ref_allele'), x.Reference_Allele = x.ref_allele; end
if isfield(x,'Tumor_Seq_Allele1') && isfield(x,'tum_allele1'), x.Tumor_Seq_Allele1 = x.tum_allele1; end
if isfield(x,'Tumor_Seq_Allele2') && isfield(x,'tum_allele2'), x.Tumor_Seq_Allele2 = x.tum_allele2; end
