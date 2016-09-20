function T = find_indel_in_GAF(GAF,F,classification,i_st,i_en,i_change,idx)
% find_indel_in_GAF(GAF,F,classification,i_st,i_en,i_change,idx)
%
% Given GAF database, annotates indels
%
%
% Alex Ramos May 2010

build = 'hg19';

if strcmpi(classification, 'Del'), is_del = true; else is_del = false; end

T.cDNAchange = repmat({'---'},length(idx),1);
T.codonchange = repmat({'---'},length(idx),1);
T.proteinchange = repmat({'---'},length(idx),1);
T.type = repmat({'---'},length(idx),1);
T.exon = nan(length(idx),1);
T.transcript = repmat({'---'},length(idx),1);
T.strand = repmat({'---'},length(idx),1);
T.gene = repmat({'---'},length(idx),1);
T.classification = repmat({'---'},length(idx),1);
T.frameness = repmat({'---'},length(idx),1);
T.idx = idx;
%T.idx = nan(length(idx),1);
T.frameness = repmat({'---'},length(idx),1);

T.refseq_mRNA_id = repmat({'---'},length(idx),1);
T.refseq_prot_id = repmat({'---'},length(idx),1);
T.swissp_acc_id = repmat({'---'},length(idx),1);
T.swissp_entry_id = repmat({'---'},length(idx),1);
T.gene_symbol = repmat({'---'},length(idx),1);
T.mRNA_id = repmat({'---'},length(idx),1);
T.description = repmat({'---'},length(idx),1);

for t=1:length(idx)
    spliceflag = false;
    init_met_stop_codon_flag = false;
    change = i_change;
    i=idx(t);
    %T.idx(t) = i;

    if isnan(i)
        T.type{t} = 'IGR';
        continue;
    end
    
    T.classification{t} = classification;
    T.change{t} = change;
    T.gene(t) = GAF.Gene(i);
    T.strand(t) = GAF.strand(i);
    T.transcript(t) = GAF.FeatureID(i);
    T.cDNApos_st(t) = NaN;
    T.cDNApos_en(t) = NaN;
    T.proteinpos_st(t) = NaN;
    T.proteinpos_en(t) = NaN;
    
	T.refseq_mRNA_id(t) = GAF.refseq_mRNA_id(i);
	T.refseq_prot_id(t) = GAF.refseq_prot_id(i);
	T.swissp_acc_id(t) = GAF.swissp_acc_id(i);
	T.swissp_entry_id(t) = GAF.swissp_entry_id(i);
	T.gene_symbol(t) = GAF.gene_symbol(i);
	T.mRNA_id(t) = GAF.mRNA_id(i);
	T.description(t) = GAF.description(i);
    
    
    plusstrand = strcmp(GAF.strand{i},'+');
    if plusstrand
        e_st = find(GAF.genomic_starts{i}-2 <= i_st & GAF.genomic_ends{i}+2 >= i_st, 1);
        e_en = find(GAF.genomic_starts{i}-2 <= i_en & GAF.genomic_ends{i}+2 >= i_en, 1);
        e_st_no_splice_sites = find(GAF.genomic_starts{i} <= i_st & GAF.genomic_ends{i} >= i_st, 1);
        e_en_no_splice_sites = find(GAF.genomic_starts{i} <= i_en & GAF.genomic_ends{i} >= i_en, 1);
    else
        e_st = find(GAF.genomic_starts{i}-2 >= i_st & GAF.genomic_ends{i}+2 <= i_st, 1);
        e_en = find(GAF.genomic_starts{i}-2 >= i_en & GAF.genomic_ends{i}+2 <= i_en, 1);
        e_st_no_splice_sites = find(GAF.genomic_starts{i} >= i_st & GAF.genomic_ends{i} <= i_st, 1);
        e_en_no_splice_sites = find(GAF.genomic_starts{i} >= i_en & GAF.genomic_ends{i} <= i_en, 1);
    end
%    e_start = find(GAF.genomic_starts{i}-2 <= R.code_start(i) & R.exon_ends{i}+2 >= R.code_start(i), 1);
%    e_end = find(GAF.genomic_starts{i}-2 <= R.code_end(i) & R.exon_ends{i}+2 >= R.code_end(i), 1);
    
    e = unique([e_st e_en]);
    %e_no_splice_sites = unique([e_st_no_splice_sites e_en_no_splice_sites]);
    if length(e) > 1, e = e(1); end

    if isempty(e)
        e = find(GAF.genomic_starts{i}-2 > i_st & GAF.genomic_ends{i}+2 < i_en, 1);
        if ~isempty(e)
            T.type{t} = 'Multi-intron_alteration';
        elseif (i_st>=GAF.tx_start(i) && i_en<=GAF.tx_stop(i)) |  (i_st<=GAF.tx_start(i) && i_en>=GAF.tx_stop(i))
            T.type{t} = 'Intron';
        elseif (i_en<=GAF.tx_start(i) && i_en<=GAF.tx_stop(i))
            if plusstrand, T.type{t} = '5''flank';
            else T.type{t} = '3''flank'; end
        elseif (i_st>=GAF.tx_start(i) && i_st>=GAF.tx_stop(i))
            if plusstrand, T.type{t} = '3''flank';
            else T.type{t} = '5''flank'; end
        else
        	error('Unable to map mutation to non-exonic site!');
        end

    else %indel falls within exon
        blank_output_flag = false;
        cpos = [i_st i_en];
        step = 1; if ~plusstrand, step = -1; end
        g_positions = [GAF.genomic_starts{i}(e):step:GAF.genomic_ends{i}(e)];
        t_positions = [GAF.transcript_starts{i}(e):GAF.transcript_ends{i}(e)];
        
        if strncmp(GAF.Gene{i},'MIR',3)     % it's a miRNA
            T.type{t} = 'miRNA';
            blank_output_flag = true;
        elseif isnan(GAF.cds_start(i))    % it's a noncoding transcript
            T.type{t} = 'Non-coding_Transcript';
            blank_output_flag = true;
        elseif isempty(e_st_no_splice_sites) | isempty(e_en_no_splice_sites)
            if (all([cpos>=GAF.tx_start(i)]) & all([cpos<=GAF.tx_stop(i)])) | (all([cpos<=GAF.tx_start(i)]) & all([cpos>=GAF.tx_stop(i)]))
                if strcmpi(classification(end-1:end), 'NP'), T.type{t} = sprintf('Splice_Site_%s',classification);
                elseif is_del, T.type{t} = 'Splice_Site_Del';
                else T.type{t} = 'Splice_Site_Ins'; end
                spliceflag = true;
				spliceinfo = get_spliceinfo(plusstrand,cpos,e,i,GAF);
				cDNA_pos = t_positions(find(g_positions == spliceinfo.cpos));
				if cDNA_pos > GAF.cds_start(i) & cDNA_pos < GAF.cds_stop(i)
					prot_pos = ceil((cDNA_pos - GAF.cds_start(i)+1)/3);
					seq = F.seq{find(strcmp(T.transcript{t},F.header))};
					prot_seq = nt2aa(seq(GAF.cds_start(i):GAF.cds_stop(i)));
					old_aa = prot_seq(prot_pos);
				end
            else
                blank_output_flag = true;
		        if (i_st<=GAF.tx_start(i) && i_st<=GAF.tx_stop(i))
		            if plusstrand, T.type{t} = '5''flank';
        		    else T.type{t} = '3''flank'; end
		        elseif (i_en>=GAF.tx_start(i) && i_en>=GAF.tx_stop(i))
        		    if plusstrand, T.type{t} = '3''flank';
	            	else T.type{t} = '5''flank'; end
    		    else
        			error('Unable to map mutation to non-exonic site!');                
	            end
	        end
        else
			cDNA_pos = [NaN NaN];
            for k=1:length(cpos)
                cDNA_pos(k) = t_positions(find(g_positions == cpos(k)));
            end
            if ~plusstrand, cDNA_pos = fliplr(cDNA_pos); end
            cDNA_change = change;
            if ~plusstrand, cDNA_change = rc(cDNA_change); end
            seq = F.seq{find(strcmp(T.transcript{t},F.header))};
            cDNA_old_base = seq(cDNA_pos(1):cDNA_pos(2));
            if is_del & ~strcmp(cDNA_old_base, cDNA_change)
                error('Error: cDNA_change does not match cDNA_old_base!');
            end
            
            if any([cDNA_pos < GAF.cds_start(i)])
                if all([cDNA_pos < GAF.cds_start(i)])  %%indel fully within UTR
                    T.type{t} = '5''-UTR';
                    blank_output_flag = true;
                    if strcmp(classification,'INS')
                        seq = [seq(1:cDNA_pos(1)) cDNA_change seq(cDNA_pos(2):end)];
                        if strfind(seq(cDNA_pos(1)-2:cDNA_pos(2)+length(change)+2),'ATG')
                            T.type{t} = 'De_novo_Start';
                        end
                    elseif strcmp(classification,'DEL')
                        if strcmp(seq(cDNA_pos(1):cDNA_pos(2)),cDNA_change)
                            seq = [seq(1:cDNA_pos(1)-1) seq(cDNA_pos(2)+1:end)];
                            if strfind(seq(cDNA_pos(1)-2:cDNA_pos(1)+2),'ATG')
                                T.type{t} = 'De_novo_Start';
                            end
                        else
                            error('Error: Genomic bases deleted does not match bases in transcript!');
                        end
                    elseif strcmpi(classification(end-1:end), 'NP')
                    	seq(cDNA_pos(1):cDNA_pos(2)) = cDNA_change;
                    	if strfind(seq(cDNA_pos(1)-2:cDNA_pos(1)+2),'ATG')
                        	T.type{t} = 'De_novo_Start';
                        end
                    end
                else
                	if strcmpi(classification(end-1:end), 'NP'), T.type{t} = sprintf('Init_Met_%s',classification);
                    elseif isdel, T.type{t} = 'Init_Met_Del';
                    else T.type{t} = 'Init_Met_Ins'; end
                    init_met_stop_codon_flag = true;
                end
            elseif any([cDNA_pos > GAF.cds_stop(i)])
                if all([cDNA_pos > GAF.cds_stop(i)])
                    T.type{t} = '3''-UTR';
                    blank_output_flag = true;
                else
                	if strcmpi(classification(end-1:end), 'NP'), T.type{t} = sprintf('Stop_Codon_%s',classification);
                    elseif is_del, T.type{t} = 'Stop_Codon_Del';
                    else T.type{t} = 'Stop_Codon_Ins'; end
                    init_met_stop_codon_flag = true;
                end
            else    
            	if strcmpi(classification(end-1:end), 'NP')
            		is_inframe = true;
                elseif mod(length([cpos(1):cpos(2)]),3) ~= 0 & is_del %frameshift deletion
                    is_inframe = false;
                    T.frameness{t} = 'out';
                elseif mod(length(change),3) ~= 0 & ~is_del    %frameshift insertion
                    is_inframe = false;
                    T.frameness{t} = 'out';
                else
                    is_inframe = true;
                    T.frameness{t} = 'in';
                end
                
                tween_codons = false;
                prot_pos = ceil((cDNA_pos - GAF.cds_start(i)+1)/3);
                prot_seq = nt2aa(seq(GAF.cds_start(i):GAF.cds_stop(i)),'ALTERNATIVESTARTCODONS',false);
                if prot_pos(1) ~= prot_pos(2)
                    tween_codons = true;  %%only matters for insertions
                    old_aa = prot_seq(prot_pos(1):prot_pos(2));
                else
                    old_aa = prot_seq(prot_pos(1));
                end
                
                [codon_start_pos, codon_end_pos] = deal(((prot_pos*3)-2)+GAF.cds_start(i)-1,(prot_pos*3)+GAF.cds_start(i)-1);
                codon_start_pos = codon_start_pos(1); codon_end_pos = codon_end_pos(2);
                old_codon = seq(codon_start_pos:codon_end_pos);                
                
                if is_inframe
                    
                    even_codon_deletion = false;
                    mut_seq = seq;
                    if strcmp(classification,'INS')
                        mut_seq = [mut_seq(1:cDNA_pos(1)) cDNA_change mut_seq(cDNA_pos(2):end)];
                        mut_prot_seq = nt2aa(mut_seq(GAF.cds_start(i):GAF.cds_stop(i)+length(change)),'ALTERNATIVESTARTCODONS',false);
                    elseif strcmp(classification,'DEL')
                        mut_seq = [mut_seq(1:cDNA_pos(1)-1) mut_seq(cDNA_pos(2)+1:end)];
                        mut_prot_seq = nt2aa(mut_seq(GAF.cds_start(i):GAF.cds_stop(i)-length(change)),'ALTERNATIVESTARTCODONS',false);
                    elseif strcmpi(classification(end-1:end), 'NP')
                    	mut_seq(cDNA_pos(1):cDNA_pos(2)) = cDNA_change;
                    	mut_prot_seq = nt2aa(mut_seq(GAF.cds_start(i):GAF.cds_stop(i)),'ALTERNATIVESTARTCODONS',false);
                    end
                    
                    if ~strcmp(old_aa, nt2aa(old_codon,'ALTERNATIVESTARTCODONS',false))
                        error('Error: old_aa variable does not match nt2aa(old_codon)');
                    end
                    
                    if is_del  %del
                        new_codon = mut_seq(codon_start_pos:codon_end_pos-length(change));
                        if isempty(new_codon)
                            even_codon_deletion = true; 
                            new_aa = '';
                        else
                            new_aa = nt2aa(new_codon,'ALTERNATIVESTARTCODONS',false);
                        end
                    elseif strcmpi(classification(end-1:end), 'NP')  %DNP
                    	new_codon = mut_seq(codon_start_pos:codon_end_pos);
                    	new_aa = nt2aa(new_codon,'ALTERNATIVESTARTCODONS',false);
                    else  %ins
                    	if tween_codons
                    		new_codon = cDNA_change;
                    	else
                    		new_codon = mut_seq(codon_start_pos:codon_end_pos+length(change));
                    	end
                        new_aa = nt2aa(new_codon,'ALTERNATIVESTARTCODONS',false);
                    end
                    
                    %%%insert checkpoint to identify if insertion is within codon or between codons
%                        mut_prot_seq = nt2aa(mut_seq(GAF.cds_start(i):GAF.cds_stop(i)));
%                        new_aa = '';
%                        for j=1:length(prot_seq)
%                            if ~strcmp(prot_seq(j),mut_prot_seq(j))
%                                if strcmp(classification,'INS')
%                                    prot_pos = [j-1 j];
%                                    new_aa = [new_aa mut_prot_seq(j)];
%                                    for k=j+1:length(mut_prot_seq)
%                                        if strcmp(mut_prot_seq(k:end), prot_seq(j:end))
%                                            break
%                                        else
%                                            new_aa = [new_aa mut_prot_seq(k)];
%                                        end
%                                    end
%                                    break
%                                elseif strcmp(classification,'DEL')
%                                
%                                end
%                            end
%                        end
%                        if isnan(prot_pos)
%                            prot_pos = ceil((cDNA_pos - GAF.cds_start(i)+1)/3);
%                        end
%                        
%                        [old_aa, new_aa] = deal(prot_seq(prot_pos), mut_prot_seq(prot_pos));
%                        [coding_seq, mut_coding_seq] = deal(seq(GAF.cds_start(i):GAF.cds_stop(i)), mut_seq(GAF.cds_start(i):GAF.cds_stop(i)));
%                        [old_codon, new_codon] = deal(coding_seq((prot_pos*3)-2:prot_pos*3), mut_coding_seq((prot_pos*3)-2:prot_pos*3));
%                        [codon_start_pos, codon_end_pos] = deal(((prot_pos*3)-2)+GAF.cds_start(i)-1,(prot_pos*3)+GAF.cds_start(i)-1);
                else
                 %%%what to do with out of frame indels?
                end
            end        
            
        end  %%end of new block
            
            
            
%        else %%cut-off old part
%            cpos = [i_st i_en];
%            if all([cpos < R.code_start(i)])
%                if plusstrand, T.type{t} = '5''-UTR'; else T.type{t} = '3''-UTR'; end
%            elseif all([cpos > R.code_end(i)])
%                if plusstrand, T.type{t} = '3''-UTR'; else T.type{t} = '5''-UTR'; end
%            else
%                T.type{t} = 'Coding';     % assume coding unless it's close to a splice-site
%                spliceflag = false;
%                init_met_stop_codon_flag = false;
%                
%                if e == e_start | e == e_end
%                    ap = cpos - R.code_start(i);
%                    bp = cpos - R.code_end(i);
%                        if any([ap >= 0 & ap <= 2]) | (ap(1) < 0 & ap(2) > 2) %indel at init. Methionine
%                        if plusstrand
%                                        iout = 'InitMet';
%                                        if is_del
%                                                T.type{t} = 'Init_Met_Del';
%                                        else
%                                                T.type{t} = 'Init_Met_Ins';
%                                        end
%                                else
%                                        iout = 'StopCodon';
%                                        if is_del
%                                                T.type{t} = 'Stop_Codon_Del';
%                                        else
%                                                T.type{t} = 'Stop_Codon_Ins';
%                                        end
%                                end
%                        cpos([cpos < R.code_start(i)]) = R.code_start(i);    %trim
%                                init_met_stop_codon_flag = true;
%                    elseif any([bp >= -2 & bp <= 0]) | (bp(1) < -2 & bp(2) > 0)
%                        if plusstrand
%                                        iout = 'StopCodon';
%                                        if is_del
%                                                T.type{t} = 'Stop_Codon_Del';
%                                        else
%                                                T.type{t} = 'Stop_Codon_Ins';
%                                        end
%                                else
%                                        iout = 'InitMet';
%                                        if is_del
%                                                T.type{t} = 'Init_Met_Del';
%                                        else
%                                                T.type{t} = 'Init_Met_Ins';
%                                        end
%                                end
%                        cpos([cpos > R.code_end(i)]) = R.code_end(i);    %trim
%                                init_met_stop_codon_flag = true;
%                    end
%                end
%        
%        %if (i_en - i_st + 1) ~= length(change)
%        %    error('"i_en - i_st + 1" should equal number of nucleotides in change variable');
%        %elseif (cpos(2)-cpos(1) + 1) ~= length(change)    % length of nucleotides in change variable differs from length of nucleotides in orf
%        %    ii = (i_st:i_en)';
%        %    change = change(find(cpos(1) == ii(:)):find(cpos(2) == ii(:)));    %trim change variable
%        %end
%        
%                if e>1     % for nonfirst exons, check left splice-site
%                    jp = cpos-R.exon_starts{i}(e);
%                    if (any([jp==-2 | jp==-1]) & all([jp < 0])) | (any([jp<0]) & any([jp>=0]))
%                        if is_del, T.type{t} = 'Splice_Site_Del'; else T.type{t} = 'Splice_Site_Ins'; end
%                        spliceflag = true;
%                        cpos = [R.exon_starts{i}(e) R.exon_starts{i}(e)];     % adjust pos to within exon for subsequent annotation
%                        if max(jp) == 0, splicedist = -1; else splicedist = max(jp); end
%                    end
%                    cpos([cpos < R.exon_starts{i}(e)]) = R.exon_starts{i}(e); %trim
%                end
%                if e<R.n_exons(i)        % for nonlast exons, check right splice-site
%                    jp = cpos-R.exon_ends{i}(e);
%                    if any([jp==2 | jp==1]) & all([jp > 0]) | (any([jp>0]) & any([jp<=0]))
%                        if is_del, T.type{t} = 'Splice_Site_Del'; else T.type{t} = 'Splice_Site_Ins'; end
%                        spliceflag = true;
%                        cpos = [R.exon_ends{i}(e) R.exon_ends{i}(e)];     % adjust pos to within exon for subsequent annotation
%                        splicedist = min(jp);
%                    end
%                    cpos([cpos > R.exon_ends{i}(e)]) = R.exon_ends{i}(e);    %trim
%                end
%
%
%                
%                T.frameness{t} = '---';
%                if mod(length([cpos(1):cpos(2)]),3) ~= 0 & is_del %frameshift deletion
%                    is_inframe = false;
%                    T.frameness{t} = 'out';
%                elseif mod(length(change),3) ~= 0 & ~is_del    %frameshift insertion
%                    is_inframe = false;
%                    T.frameness{t} = 'out';
%                else
%                    is_inframe = true;
%                    T.frameness{t} = 'in';
%                end
%                % For coding and splice-site mutations: find position in cDNA and protein
%                
%                orf = [];
%                coding_exon = 0;
%                nframeshifts = 0;
%
%                if plusstrand, forfrom=1; forstep=+1; forto=R.n_exons(i);
%                else forfrom=R.n_exons(i); forstep=-1; forto=1;
%                end
%
%                for e=forfrom:forstep:forto
%                    st = R.exon_starts{i}(e);
%                    en = R.exon_ends{i}(e);
%                    fr = R.exon_frames{i}(e);
%                    
%                    if ~isempty(orf) && fr~=-1 && mod(length(orf),3)~=fr
%                        % recover from known programmed frameshift somewhere in the preceding exon
%                        nframeshifts = nframeshifts+1;
%                        pfs = mod(fr-length(orf),3);
%                        if pfs==1, orf(end+1) = 'A';
%                        else orf(end) = [];
%                        end
%                    end
%                    
%            %%%%%%%%%%%%!A - builds orf and identifies orf position of mutation
%                    if st<=R.code_end(i) && en>=R.code_start(i)
%                        coding_exon = coding_exon + 1;
%                        if R.code_start(i)>st, st = R.code_start(i); end
%                        if R.code_end(i)<en, en = R.code_end(i); end
%                        if all([cpos >= st & cpos <= en])    % indel is within this exon
%                            T.exon(t) = coding_exon;
%                            if plusstrand
%                orf_pos = length(orf) + (cpos-st+1);
%                            else 
%                orf_pos = length(orf) + (en-cpos+1);
%                orf_pos = fliplr(orf_pos);
%                            end
%                            orf_codon_start_pos = orf_pos - mod(orf_pos-1,3);
%                            orf_codon_end_pos = orf_codon_start_pos + 2;
%                        end
%                        d = upper(genome_region(R.chr(i),st,en,build));
%                        if ~plusstrand, d = rc(d); end
%                        orf = [orf d];
%                    end
%
%                end
%                
%                T.cDNApos_st(t) = orf_pos(1);
%                T.cDNApos_en(t) = orf_pos(2);
%                aa_number_st = orf_codon_end_pos(1)/3;
%                aa_number_en = orf_codon_end_pos(2)/3;
%                T.proteinpos_st(t) = aa_number_st;
%                T.proteinpos_en(t) = aa_number_en;
%
%                if orf_codon_end_pos(2) > length(orf)     % problem!    truncated ORF
%                    orf = [orf 'AA'];
%                end
%
%%                if verbosity>=2
%%                    T.orf{t} = orf;
%%                    T.protein{t} = my_nt2aa(orf);
%%                    T.nframeshifts(t) = nframeshifts;
%%                end
%        
%                old_codon = orf(orf_codon_start_pos(1):orf_codon_end_pos(2));
%                old_aa = my_nt2aa(old_codon);
%                if ~spliceflag
%                    if is_inframe
%                        even_codon_deletion = false;    
%                        if is_del
%                            c_del = orf(orf_pos(1):orf_pos(2));
%                            %if ~plusstrand, c_del = rc(c_del); end
%                            new_codon = regexprep(old_codon, c_del, '');
%                            if strcmp(new_codon, ''), even_codon_deletion = true; end
%                        elseif ~is_del    %insertion
%                            if orf_codon_start_pos(1) == orf_codon_start_pos(2)    %insertion falls within a codon
%                                if ~plusstrand, change1 = rc(change); else change1 = change; end        
%                                if orf_pos(1) - orf_codon_start_pos(1) == 1    %insertion falls within 2nd&3rd bases of codon
%                                    new_codon = [old_codon(1:2) change1 old_codon(end)];
%                                elseif orf_pos(1) - orf_codon_start_pos(1) == 0    %insertion falls within 1st&2nd bases of codon
%                                    new_codon = [old_codon(1) change1 old_codon(2:end)];
%                                end
%                                tween_codons = false;
%                            else    %insertion falls between two codons
%                                tween_codons = true;
%                                if plusstrand, new_codon = change; else new_codon = rc(change); end
%                                %next_aa = my_nt2aa(orf(orf_codon_start_pos(1)+3:orf_codon_end_pos(2)+3));
%                                next_aa = old_aa(end);
%                                next_aa_pos = aa_number_st + 1;
%                             end
%                        end
%                        new_aa = my_nt2aa(new_codon);
%                    end
%                end
                
                
                    
        if spliceflag
            T.cDNAchange{t} = sprintf('c.%d_splice',cDNA_pos);
            %if ~plusstrand, splicedist = -splicedist; end
            T.codonchange{t} = sprintf('c.e%d%+d',e,spliceinfo.splicedist);
            if exist('old_aa','var')
            	T.proteinchange{t} = sprintf('p.%s%d_splice',old_aa,prot_pos);
            end
        elseif init_met_stop_codon_flag
            T.cDNAchange{t} = sprintf('c.%d_%s',cDNA_pos(1),iout);
            %T.codonchange{t} = sprintf('c.e%d_%s',e,iout);
            %T.proteinchange{t} = sprintf('p.%s%d_%s',old_aa,aa_number_st,iout);
        elseif blank_output_flag
        	%pass
        elseif ~is_inframe
            if strcmpi(classification, 'ins')
                T.type{t} = 'Frame_Shift_Ins';
                %if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%dins%s',orf_pos(1),orf_pos(2),rc(change));
                %else T.cDNAchange{t} = sprintf('c.%d_%dins%s',orf_pos(1),orf_pos(2),change);
                %end
                T.cDNAchange{t} = sprintf('c.%d_%dins%s',cDNA_pos(1),cDNA_pos(2),cDNA_change);
            elseif strcmpi(classification, 'del')
                T.type{t} = 'Frame_Shift_Del';
                %if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),rc(change));
                %else T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),change);
                %end
                %if ~plusstrand cc = rc(change);, else cc = change; end
                if cDNA_pos(1) == cDNA_pos(2)
                    T.cDNAchange{t} = sprintf('c.%ddel%s',cDNA_pos(1),cDNA_change);
                else
                    T.cDNAchange{t} = sprintf('c.%d_%ddel%s',cDNA_pos(1),cDNA_pos(2),cDNA_change);
                end
        
            end
            T.codonchange{t} = sprintf('c.(%d-%d)%sfs',codon_start_pos,codon_end_pos,old_codon);
            T.proteinchange{t} = sprintf('p.%s%dfs',old_aa(1),prot_pos(1));
        elseif is_inframe % Coding
        	if strcmpi(classification(end-1:end), 'NP')

	            if old_aa==new_aa, T.type{t} = 'Synonymous';
	            elseif new_aa=='*', T.type{t} = 'Nonsense';
	            elseif old_aa=='*', T.type{t} = 'Read-through';
	            else T.type{t} = 'Missense';
	            end

		        T.cDNAchange{t} = sprintf('c.%d_%d%s>%s',cDNA_pos(1),cDNA_pos(2),cDNA_old_base,cDNA_change);
	            T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',codon_start_pos,codon_end_pos,old_codon,new_codon);
	            if prot_pos(1) == prot_pos(2)
		            T.proteinchange{t} = sprintf('p.%s%d%s',old_aa,prot_pos(1),new_aa);
		        else
		        	T.proteinchange{t} = sprintf('p.%d_%d%s>%s',prot_pos(1),prot_pos(2),old_aa,new_aa);
		        end

            elseif ~is_del
                T.type{t} = 'In_Frame_Ins';
                %if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%dins%s' ,orf_pos(1),orf_pos(2),rc(change));
                %else T.cDNAchange{t} = sprintf('c.%d_%dins%s',orf_pos(1),orf_pos(2),change);
                %end
        		T.cDNAchange{t} = sprintf('c.%d_%dins%s',cDNA_pos(1),cDNA_pos(2),cDNA_change);
                if tween_codons
                    T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',codon_start_pos,codon_end_pos,old_codon,[old_codon(1:3) new_codon old_codon(4:end)]);
                    T.proteinchange{t} = sprintf('p.%s%d_%s%dins%s',old_aa(1),prot_pos(1),old_aa(2),prot_pos(2),new_aa);
                else
                    T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',codon_start_pos,codon_end_pos,old_codon,new_codon);
                    T.proteinchange{t} = sprintf('p.%s%d>%s',old_aa,prot_pos(1),new_aa);
                end
            elseif is_del
                T.type{t} = 'In_Frame_Del';
                %if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),rc(c_del));
                %else T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),c_del);
                %end
                if cDNA_pos(1) == cDNA_pos(2)
                    T.cDNAchange{t} = sprintf('c.%ddel%s',cDNA_pos(1),cDNA_old_base);
                else
                    T.cDNAchange{t} = sprintf('c.%d_%ddel%s',cDNA_pos(1),cDNA_pos(2),cDNA_old_base);
                end
        
                if even_codon_deletion
                    T.codonchange{t} = sprintf('c.(%d-%d)%sdel',codon_start_pos,codon_end_pos,old_codon);
                    if length(old_aa) > 1
                        T.proteinchange{t} = sprintf('p.%s%d_%s%ddel',old_aa(1),prot_pos(1),old_aa(end),prot_pos(2));
                    else
                        T.proteinchange{t} = sprintf('p.%s%ddel',old_aa,prot_pos(1));
                    end
                else
                    T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',codon_start_pos,codon_end_pos,old_codon,new_codon);
                    if length(old_aa) > 1
                        T.proteinchange{t} = sprintf('p.%s%d_%s%d>%s',old_aa(1),prot_pos(1),old_aa(end),prot_pos(2),new_aa);
                    else
                        T.proteinchange{t} = sprintf('p.%s%d>%s',old_aa,prot_pos(1),new_aa);
                    end
                end        
            end
        end
%            end
%        end
%    end
    

end     % next transcript
end  % end of main function

function spliceinfo = get_spliceinfo(plusstrand,cpos,e,i,GAF)
	spliceinfo = [];
	if plusstrand
		if any(cpos < GAF.genomic_starts{i}(e)) %% indel on left side of exon
			spliceinfo.splicedist = min(cpos(find(cpos < GAF.genomic_starts{i}(e))) - GAF.genomic_starts{i}(e));
			spliceinfo.cpos = GAF.genomic_starts{i}(e);
		elseif any(cpos > GAF.genomic_ends{i}(e)) %% indel on right side of exon
			spliceinfo.splicedist = min(cpos(find(cpos > GAF.genomic_ends{i}(e))) - GAF.genomic_ends{i}(e));
			spliceinfo.cpos = GAF.genomic_ends{i}(e);
		end
	else
		if any(cpos > GAF.genomic_starts{i}(e)) %% indel on right side of exon (on + strand)
			spliceinfo.splicedist = min(cpos(find(cpos > GAF.genomic_starts{i}(e))) - GAF.genomic_starts{i}(e));
			spliceinfo.cpos = GAF.genomic_starts{i}(e);
		elseif any(cpos < GAF.genomic_ends{i}(e)) %% indel on left side of exon (on + strand)
			spliceinfo.splicedist = min(cpos(find(cpos < GAF.genomic_ends{i}(e))) - GAF.genomic_ends{i}(e));
			spliceinfo.cpos = GAF.genomic_ends{i}(e);
		end
		spliceinfo.splicedist = -spliceinfo.splicedist;
	end

