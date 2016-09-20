function M = classify_muts_GAF(M1, P, GAF, F)
%
% given a struct with fields:
%     chr, start, end, change
%
% finds the highest-priority change caused by the specified change
%     "change" must have one of the following formats:
%                A         point mutation
%                C/G     point mutation to two different alleles
%                -A        deletion
%                +CC     insertion
%
% and returns a struct with the added fields:
%     gene, type, ridx(index in RefSeq), transcript, strand,
%     genomechange, cDNAchange, codonchange, proteinchange
%


%%%%%%%%%%%%%%%!@ Add new option to output exon coords
if ~exist('P','var'), P=[]; end
%P = impose_default_value(P,'build','hg18');
%P = impose_default_value(P,'build_refseq',P.build);
%P = impose_default_value(P,'build_genome_region',P.build);
P = impose_default_value(P,'impute_promoters',true);
P = impose_default_value(P,'imputed_promoter_size',3000);
%P = impose_default_value(P,'include_new_fields',false);
%P = impose_default_value(P,'include_sequences',false);
P = impose_default_value(P,'include_all_transcripts',false);
%P = impose_default_value(P,'include_exon_coordinates',false);

%verbosity = 0;
%if P.include_new_fields, verbosity = 1; end
%if P.include_sequences, verbosity = 2; end

require_fields(M1,{'chr','start','end','change'});
nm = slength(M1);

type_priority = {'Non-mutation','De_novo_Start','Nonsense','Read-through','Missense',...
    'Splice_site','Synonymous','Promoter','miRNA','miRNA_vicinity','3''-UTR','5''-UTR',...
    'Intron','5''flank','3''flank','Non-coding_Transcript','IGR'};
    
indel_type_priority = {'De_novo_Start','In_Frame_Del','In_Frame_Ins','Frame_Shift_Del','Frame_Shift_Ins',...
    'Init_Met_Del','Init_Met_Ins','Stop_Codon_Del','Stop_Codon_Ins','Stop_codon_indel',...
    'Splice_Site_Del','Splice_Site_Ins','5''-UTR','3''-UTR','Promoter','miRNA','Intron',...
    '5''flank','3''flank','Non-coding_Transcript','IGR'};

dnp_type_priority = {'De_novo_Start','Nonsense','Read-through','Missense','Init_Met_DNP','Init_Met_TNP',...
    'Init_Met_ONP','Stop_Codon_DNP','Stop_Codon_TNP','Stop_Codon_ONP','Splice_Site_DNP',...
    'Splice_Site_TNP','Splice_Site_ONP','Synonymous','miRNA','3''-UTR','5''-UTR','Promoter',...
    'Intron','5''flank','3''flank','Non-coding_Transcript','IGR'};

kflds = {'idx','gene','strand','classification','type','genomechange','transcript',...
    'cDNAchange','codonchange','proteinchange','other_transcripts',...
    'refseq_mRNA_id','refseq_prot_id','swissp_acc_id','swissp_entry_id','gene_symbol',...
    'mRNA_id','description'};
%
%if P.include_new_fields
%    kflds = [kflds 'exon','cDNApos','proteinpos'];
%end
%if P.include_sequences
%    kflds = [kflds 'nframeshifts' 'tx_seq' 'orf' 'protein'];
%end
%if P.include_exon_coordinates
%    kflds = [kflds 'n_exons' 'exon_coords'];
%end

if nm==0
    fprintf('Adding blank fields to blank struct.\n');
    M = M1;
    for i=1:length(kflds), M = setfield(M,kflds{i},[]); end
    return
end

% load databases of genomic features
%%%%%%%%%%%!@ Change this part to use new UCSC file
%if ~exist('R','var')
%    R = load_refseq(P.build_refseq);
%    fprintf('Forcing TP53 annotations to transcript NM_000546 (full-length isoform)\n');
%    i1 = grep('^TP53$',R.gene,1);
%    i2 = grep('NM_000546',R.transcript,1);
%    i3 = setdiff(i1,i2);
%    keep = setdiff(1:slength(R),i3);
%    R = reorder_struct(R,keep);
%end

% ADD PROMOTERS TO DATABASE REGIONS
GAF.gene_start = GAF.tx_start;
GAF.gene_end = GAF.tx_stop;
if P.impute_promoters
    idx = find(strcmp('+',GAF.strand));
    GAF.gene_start(idx) = GAF.gene_start(idx) - P.imputed_promoter_size;
    idx = find(strcmp('-',GAF.strand));
    GAF.gene_start(idx) = GAF.gene_start(idx) + P.imputed_promoter_size;
end

% CLASSIFY EACH MUTATION

M=[];

if nm>10000, step=1000; else step=100; end
for i=1:nm
    if ~mod(i,step), fprintf('%d/%d ',i,nm); end

    change = M1.change{i};
    if length(change)==3 && change(2)=='/'     % point mutation to two different alleles: try both
        changes = {change(1);change(3)};
    else
        changes = {change};
    end

    T = cell(length(changes),1);
    for x=1:length(changes)
        T{x} = find_mut_in_GAF(GAF,F,M1.chr{i},M1.start(i),M1.end(i),changes{x},[]);
    end
    T = concat_structs(T);

    if P.include_all_transcripts
        idx = 1:slength(T);
    elseif strncmp(change, '+', 1) | strncmp(change, '-', 1)    %indel
        idx = filter_idxs(T, GAF, indel_type_priority);
    elseif strncmp(change, '~', 1)
        idx = filter_idxs(T, GAF, dnp_type_priority);
    else
        idx = filter_idxs(T, GAF, type_priority);
        
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
    end
    if iscell(o_string)
        o_string = o_string{1};
        o_string = o_string(1:end-1);
    end
    
    N2.other_transcripts = {o_string};
    
    %%%%%%%%%%%%
    n_exons = GAF.n_exons(idx);    


    N2.n_exons = n_exons;
    %%%%%%%%%%%%%%%!@
    
    
    N2 = keep_fields(N2,kflds);
    if slength(N1) ~= slength(N2)    %%temporary fix for super-rare N2s that have an extra index value
        N2.idx = N2.idx(1:end-1);
    end
    
    M{i} = merge_structs({N1,N2});

end     % next mutation

M = concat_structs(M);
end

function idx = filter_idxs(T, GAF, priority_list)
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
    if ~all(isnan(GAF.code_len(filtered_idxs)))
      %%Only consider transcripts with longest coding length
      ii = find(GAF.code_len(filtered_idxs) == max(GAF.code_len(filtered_idxs)));
      filtered_idxs = filtered_idxs(ii);

    end
    
    filtered_idxs = filtered_idxs(1);
    idx = find(T.idx == filtered_idxs);

end
