function M = classify_muts(M1, P, R)
%
% given a struct with fields:
%   chr, start, end, change
%
% finds the highest-priority change caused by the specified change
%   "change" must have one of the following formats:
%        A     point mutation
%        C/G   point mutation to two different alleles
%        -A    deletion
%        +CC   insertion
%
% and returns a struct with the added fields:
%   gene, type, ridx(index in RefSeq), transcript, strand,
%   genomechange, cDNAchange, codonchange, proteinchange
%


%%%%%%%%%%%%%%%!@ Add new option to output exon coords
if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'build_refseq',P.build);
P = impose_default_value(P,'build_genome_region',P.build);
P = impose_default_value(P,'impute_promoters',true);
P = impose_default_value(P,'imputed_promoter_size',3000);
P = impose_default_value(P,'include_new_fields',false);
P = impose_default_value(P,'include_sequences',false);
P = impose_default_value(P,'include_all_transcripts',false);
P = impose_default_value(P,'include_exon_coordinates',false);
P = impose_default_value(P,'force_manual_transcript_choices',false);

verbosity = 0;
if P.include_new_fields, verbosity = 1; end
if P.include_sequences, verbosity = 2; end

require_fields(M1,{'chr','start','end','change'});
nm = slength(M1);

type_priority = {'Non-mutation','De_novo_Start','Nonsense','Read-through','Missense','Splice_site',...
  'Synonymous','Promoter','miRNA','miRNA_vicinity','3''-UTR','5''-UTR','Intron','Non-coding_Transcript','IGR'};
  
indel_type_priority = {'In_Frame_Del','In_Frame_Ins','Frame_Shift_Del','Frame_Shift_Ins',...
	'Init_Met_Del','Init_Met_Ins','Stop_Codon_Del','Stop_Codon_Ins','Stop_codon_indel','Splice_Site_Del','Splice_Site_Ins','5''-UTR','3''-UTR','Promoter','miRNA','Intron','Non-coding_Transcript','IGR'};

dnp_type_priority = {'Nonsense','Read-through','Missense','Init_Met_DNP','Init_Met_TNP','Init_Met_ONP','Stop_Codon_DNP',...
	'Stop_Codon_TNP','Stop_Codon_ONP','Splice_Site_DNP','Splice_Site_TNP','Splice_Site_ONP','Synonymous','miRNA','3''-UTR','5''-UTR','Promoter','Intron','Non-coding_Transcript','IGR'};

kflds = {'idx','gene','strand','classification','type','genomechange','transcript','cDNAchange','codonchange', ...
         'proteinchange','other_transcripts'};

if P.include_new_fields
  kflds = [kflds 'exon','cDNApos','proteinpos'];
end
if P.include_sequences
  kflds = [kflds 'nframeshifts' 'tx_seq' 'orf' 'protein'];
end
if P.include_exon_coordinates
  kflds = [kflds 'n_exons' 'exon_coords'];
end

if nm==0
  fprintf('Adding blank fields to blank struct.\n');
  M = M1;
  for i=1:length(kflds), M = setfield(M,kflds{i},[]); end
  return
end

% load databases of genomic features
%%%%%%%%%%%!@ Change this part to use new UCSC file
if ~exist('R','var')
  R = load_refseq(P.build_refseq);
end

if P.force_manual_transcript_choices
  fprintf('Forcing TP53 annotations to transcript NM_000546 (full-length isoform)\n');
  i1 = grep('^TP53$',R.gene,1);
  i2 = grep('NM_000546',R.transcript,1);
  i3 = setdiff(i1,i2);
  keep = setdiff(1:slength(R),i3);
  R = reorder_struct(R,keep);

  fprintf('Forcing TTN annotations to transcript NM_133379\n');
  i1 = grep('^TTN$',R.gene,1);
  i2 = grep('NM_133379',R.transcript,1);
  i3 = setdiff(i1,i2);
  keep = setdiff(1:slength(R),i3);
  R = reorder_struct(R,keep);
end

% ADD PROMOTERS TO DATABASE REGIONS
R.gene_start = R.tx_start;
R.gene_end = R.tx_end;
if P.impute_promoters
  idx = find(strcmp('+',R.strand));
  R.gene_start(idx) = R.gene_start(idx) - P.imputed_promoter_size;
  idx = find(strcmp('-',R.strand));
  R.gene_end(idx) = R.gene_end(idx) + P.imputed_promoter_size;
end

% CLASSIFY EACH MUTATION

M=[];

if nm>10000, step=1000; else step=100; end
for i=1:nm
  if ~mod(i,step), fprintf('%d/%d ',i,nm); end

  change = M1.change{i};
  if length(change)==3 && change(2)=='/'   % point mutation to two different alleles: try both
    changes = {change(1);change(3)};
  else
    changes = {change};
  end

  T = cell(length(changes),1);
  for x=1:length(changes)
    if isnumeric(M1.chr)
      T{x} = find_mut_in_refseq(R,P.build_genome_region,M1.chr(i),M1.start(i),M1.end(i),changes{x},[],verbosity);
    else
      T{x} = find_mut_in_refseq(R,P.build_genome_region,M1.chr{i},M1.start(i),M1.end(i),changes{x},[],verbosity);
    end
  end
  T = concat_structs(T);
  
  %if i == 5
  %	error('test');
  %end

%%%Checks for presence of NR_xxx (non-coding) transcripts and removes them only if an NM_xxx transcript is in T as well
%  if ~P.include_all_transcripts	
%	if slength(T) > 1
%	  if any(strncmp(T.transcript, 'NR_', 3)) & any(strncmp(T.transcript, 'NM_', 3))
%		NM_idx = strncmp(T.transcript, 'NM_', 3);
%		T = reorder_struct(T,NM_idx);
%	  end
%	end
%  end

  if P.include_all_transcripts
    idx = 1:slength(T);
  elseif strncmp(change, '+', 1) | strncmp(change, '-', 1)  %indel
    idx = filter_idxs(T, R, indel_type_priority);
%  	for p=1:length(indel_type_priority)
%  		idx = find(strcmpi(T.type,indel_type_priority{p}),1);
%  		if ~isempty(idx), break; end
%  	end
  elseif strncmp(change, '~', 1)
    idx = filter_idxs(T, R, dnp_type_priority);
%  	for p=1:length(dnp_type_priority)
%  		idx = find(strcmp(T.type,dnp_type_priority{p}),1);
%  		if ~isempty(idx), break; end
%  	end
  else
    idx = filter_idxs(T, R, type_priority);
%    for p=1:length(type_priority)
%      idx = find(strcmpi(T.type,type_priority{p}));
%      if ~isempty(idx), break; end
%    end
%    
%    %%Only consider transcripts with longest coding length
%    idx = find(R.code_len(T.idx(idx)) == max(R.code_len(T.idx(idx))));
%    
%    %%Only consider transcript with high version order
%    nn = zeros(length(idx),1);
%    for j=1:length(idx)
%      nn(j) = str2num(T.transcript{idx(j)}(end));
%    end
%    [nn, ni] = sort(nn,'descend');
%    idx = idx(ni);
%    idx = idx(1);
    
  end

  if isempty(idx), error('idx should not be empty'); end

  N1 = reorder_struct(M1,i*ones(length(idx),1));
  N2 = reorder_struct(T,idx);
  
  %%%%%%%add string containing other transcripts and protein changes
  o_string = '';
  b_idx = setdiff(1:slength(T), idx);
  for z=b_idx
    if strcmp(T.proteinchange(z),'---')
    	o_string = strcat(o_string,strcat(T.transcript(z),'_',T.type(z)),'|');
    else
       	o_string = strcat(o_string,strcat(T.transcript(z),'_',T.type(z),'_',T.proteinchange(z)),'|');
    end
  	% o_string = strcat(o_string,strcat(T.transcript(z),'_',T.proteinchange(z)),'|');
  end
  if iscell(o_string)
  	o_string = o_string{1};
  	o_string = o_string(1:end-1);
  end
  
  N2.other_transcripts = {o_string};
  
  %%%%%%%%%%%%
  n_exons = R.n_exons(idx);	
  if P.include_exon_coordinates
  	idxe = zeros(slength(N2),1);
  	for j=1:slength(N2)
		idxe(j) = find(strcmp(T.transcript(j), R.transcript),1);
	end


	zz = [R.exon_starts(idxe) R.exon_ends(idxe)];

	coords_cat = cell(size(zz,1),1);
	for k=1:size(zz,1)
		coords = cellstr(strcat(num2str(zz{k,1}), repmat('-',length(zz{k,1}),1), num2str(zz{k,2})));
		ll = '';
		for q=1:length(coords)
			if q == 1
				ll = strcat(ll, coords(q));
			else
				ll = strcat(ll, '|', coords(q));
			end
		end
		coords_cat(k) = ll;		
	end
	N2.exon_coords = coords_cat;
  end

  N2.n_exons = n_exons;
  %%%%%%%%%%%%%%%!@
  
  
  N2 = keep_fields(N2,kflds);
  if slength(N1) ~= slength(N2)  %%temporary fix for super-rare N2s that have an extra index value
  	N2.idx = N2.idx(1:end-1);
  end
  M{i} = merge_structs({N1,N2});

end   % next mutation

M = concat_structs(M);
end

%function idx = filter_idxs(T, R, priority_list)
%	filtered_idxs = [];
%  	idx = [];
%  	%ii=
%    %%Only consider transcripts with highest position on priority list
%    for p=1:length(priority_list)
%      idx = find(strcmpi(T.type,priority_list{p}));
%      if ~isempty(idx), break; end
%    end
%    
%    
%    %Skip IGR muts
%    if ~isnan(T.idx) & ~all(isnan(R.code_len(filtered_idxs)))
%  	  filtered_idxs = T.idx(idx);
%      %%Only consider transcripts with longest coding length
%      ii = find(R.code_len(filtered_idxs) == max(R.code_len(filtered_idxs)));
%      filtered_idxs = filtered_idxs(ii);
%      
%      filtered_idxs = filtered_idxs(1);
%      idx = find(T.idx == filtered_idxs);
%      %%Only consider transcript with high version order
%      %nn = zeros(length(idx),1);
%      %for j=1:length(ii)
%      %  nn(j) = str2num(T.transcript{idx(j)}(end));
%      %end
%      %[nn, ni] = sort(nn,'descend');
%      %idx = idx(ni);
%    end
%        
%    %idx = idx(1);
%    
%    if length(idx) > 1
%    	idx = idx(1);
%    end 
%    
%end

function idx = filter_idxs(T, R, priority_list)
	if all(strcmp(T.gene,'---')) | all(strcmp(T.gene,'---'))
		idx = 1;
		return
	elseif all(strcmp(T.gene,''))
		filtered_idxs = T.idx;
	else
		filtered_idxs = T.idx(~strcmp(T.gene,''));  %filter out transcripts with no gene
	end
  	idx = [];
  	%ii=
    %%Only consider transcripts with highest position on priority list
    types = cell(length(filtered_idxs),1);
    for f=1:length(filtered_idxs)
    	types(f) = T.type(find(T.idx == filtered_idxs(f)));
    end
    for p=1:length(priority_list)
      idx = find(strcmpi(types, priority_list{p}));
      if ~isempty(idx), break; end
    end
    filtered_idxs = filtered_idxs(idx); 
    
    %Skip IGR muts and noncoding transcripts
    if ~all(isnan(R.code_len(filtered_idxs)))
      %filtered_idxs = T.idx(idx);
      %%Only consider transcripts with longest coding length
      ii = find(R.code_len(filtered_idxs) == max(R.code_len(filtered_idxs)));
      filtered_idxs = filtered_idxs(ii);
      

      
      %%Only consider transcript with high version order
      %nn = zeros(length(idx),1);
      %for j=1:length(ii)
      %  nn(j) = str2num(T.transcript{idx(j)}(end));
      %end
      %[nn, ni] = sort(nn,'descend');
      %idx = idx(ni);
    end
    
    filtered_idxs = filtered_idxs(1);
    idx = find(T.idx == filtered_idxs);
    %idx = idx(1);
    %if length(idx) > 1
    %	idx = idx(1);
    %end 
end
