function T = find_mut_in_GAF(GAF,F,chr,st,en,change,idx)
% find_mut_in_refseq(R, build, chromosome, start, end, change, idx, verbosity)
%
% "change" must have one of the following formats:
%                                A                 point mutation
%                                -A                deletion
%                                +CC         insertion
%                                ~TT         DNP
%
% given RefSeq database in struct R (loaded using load_refseq),
% finds all Refseq transcripts overlapping the indicated mutation,
% and determines the outcome of mutation to the nucleotide specified in "change"
%
% can force annotation to a particular transcript using "idx"
%
% (Mike Lawrence 2008-2010)
% Alex Ramos 


if ~isfield(GAF,'gene_start'), GAF.gene_start = GAF.tx_start; end
if ~isfield(GAF,'gene_end'), GAF.gene_end = GAF.tx_stop; end

if isnumeric(GAF.chr) && isnumeric(chr)
        % OK
elseif ~isnumeric(GAF.chr) && ~isnumeric(chr) && (strncmp(chr,'chr',3)==strncmp(GAF.chr{1},'chr',3))
        % OK
else         % try to reconcile incompatible chr formats
        GAF.chr = convert_chr(GAF.chr);
        chr = convert_chr(chr);
end

pos = st;        %variable for point mutations
is_nonSNP = false;
build = 'hg19';
old_base = upper(genome_region(chr,st,en,build));
change = upper(change);
non_mut_flag = strcmp(old_base,change);

%if ~exist('verbosity','var'), verbosity = 0; end         % flag whether to return additional info

if length(change) > 1 & (strncmp(change,'+',1) | strncmp(change,'-',1) | strncmp(change,'~',1)) %is indel or DNP
    is_nonSNP = true;
    %O = struct('adjust_coordinates_for_Broad',false);        %for TCGA ovarian indel testing
    I = preprocess_nonSNP_call_GAF(GAF,chr,st,en,change);
    change = I.change;
    idx = I.idx;
    
elseif length(change)==1 | length(change) == 2  
        if ~exist('idx','var') || isempty(idx)
                if isnumeric(GAF.chr)
                        chridx = find(GAF.chr==chr);
                else
                        chridx = find(strcmp(GAF.chr,chr));
                end
                %if isempty(chridx), fprintf('Warning: unable to reconcile chromosome formats\n'); end
                idx_plus_strand = chridx(GAF.gene_start(chridx)<=pos & GAF.gene_end(chridx)>=pos);
                idx_minus_strand = chridx(GAF.gene_start(chridx)>=pos & GAF.gene_end(chridx)<=pos);
                idx = [idx_plus_strand; idx_minus_strand];
        end
else
                 error('Invalid "change" = %s',change);
end

if isempty(idx)         % IGR mutation
        idx = nan;
end

T=[];
T.idx = idx;
T.gene = repmat({'---'},length(idx),1);
T.strand = repmat({'---'},length(idx),1);
T.transcript = repmat({'---'},length(idx),1);
T.type = repmat({'---'},length(idx),1);
T.exon = nan(length(idx),1);
T.cDNApos = nan(length(idx),1);
T.proteinpos = nan(length(idx),1);

if is_nonSNP
    T.classification = repmat({I.classification},length(idx),1);
    T.genomechange = repmat({I.genomechange},length(idx),1);
else
    T.classification = repmat({'SNP'},length(idx),1);
    
    if iscell(chr)
        error('chr should not be a cell');
    elseif isnumeric(chr)
        l = chrlabel(chr);
        T.genomechange = repmat({['g.' l{1} ':'...
        num2str(pos) old_base '>' change]},length(idx),1);
    elseif strncmp(chr,'chr',3)
        T.genomechange = repmat({['g.' chr ':'...
        num2str(pos) old_base '>' change]},length(idx),1);
    else
        T.genomechange = repmat({['g.chr' chr ':'...
        num2str(pos) old_base '>' change]},length(idx),1);
    end
end

T.cDNAchange = repmat({'---'},length(idx),1);
T.codonchange = repmat({'---'},length(idx),1);
T.proteinchange = repmat({'---'},length(idx),1);

T.refseq_mRNA_id = repmat({'---'},length(idx),1);
T.refseq_prot_id = repmat({'---'},length(idx),1);
T.swissp_acc_id = repmat({'---'},length(idx),1);
T.swissp_entry_id = repmat({'---'},length(idx),1);
T.gene_symbol = repmat({'---'},length(idx),1);
T.mRNA_id = repmat({'---'},length(idx),1);
T.description = repmat({'---'},length(idx),1);

%if verbosity>=2
%        T.tx_seq = repmat({'---'},length(idx),1);
%        T.orf = repmat({'---'},length(idx),1);
%        T.protein = repmat({'---'},length(idx),1);
%        T.nframeshifts = nan(length(idx),1);
%end
spliceflag = false;
%if st == 66455993
%	error('LOL')
%end


if is_nonSNP
    %if strcmpi(I.classification(end-1:end), 'NP')        %DNP,TNP, or ONP
    %    I2 = find_DNP_in_refseq(R,build,I.classification,chr,I.pos1,I.pos2,change,idx,verbosity);
    %elseif strcmpi(I.classification, 'DEL') | strcmpi(I.classification, 'INS')
        I2 = find_indel_in_GAF(GAF,F,I.classification,I.pos1,I.pos2,change,idx);
    	
    %else
    %    error('Multi-nucleotide variant does not match valid types, should not happen!');
    %end
    T.gene = I2.gene;
    T.strand = I2.strand;
    T.transcript = I2.transcript;
    T.type = I2.type;
    T.exon = I2.exon;
    T.cDNAchange = I2.cDNAchange;
    T.codonchange = I2.codonchange;
    T.proteinchange = I2.proteinchange;
    
	T.refseq_mRNA_id = I2.refseq_mRNA_id;
	T.refseq_prot_id = I2.refseq_prot_id;
	T.swissp_acc_id = I2.swissp_acc_id;
	T.swissp_entry_id = I2.swissp_entry_id;
	T.gene_symbol = I2.gene_symbol;
	T.mRNA_id = I2.mRNA_id;
	T.description = I2.description;
else

    for t=1:length(idx)

        i=idx(t);
        spiceflag = false;

        if isnan(i)
            T.type{t} = 'IGR';
            continue;
        elseif non_mut_flag
   			T.type{t} = 'Non-mutation';
   			continue;
        end

        T.gene(t) = GAF.Gene(i);
        T.strand(t) = GAF.strand(i);
        T.transcript(t) = GAF.FeatureID(i);
        
		T.refseq_mRNA_id(t) = GAF.refseq_mRNA_id(i);
		T.refseq_prot_id(t) = GAF.refseq_prot_id(i);
		T.swissp_acc_id(t) = GAF.swissp_acc_id(i);
		T.swissp_entry_id(t) = GAF.swissp_entry_id(i);
		T.gene_symbol(t) = GAF.gene_symbol(i);
		T.mRNA_id(t) = GAF.mRNA_id(i);
		T.description(t) = GAF.description(i);

        plusstrand = strcmp(GAF.strand{i},'+');

        %if verbosity>=2, T.tx_seq{t} = sub_get_tx_seq(i); end

        %!A - indel
        if plusstrand
            e = find(GAF.genomic_starts{i}-2 <= pos & GAF.genomic_ends{i}+2 >= pos, 1);
            e_no_splice_sites = find(GAF.genomic_starts{i} <= pos & GAF.genomic_ends{i} >= pos, 1);
        else
            e = find(GAF.genomic_starts{i}+2 >= pos & GAF.genomic_ends{i}-2 <= pos, 1);
            e_no_splice_sites = find(GAF.genomic_starts{i} >= pos & GAF.genomic_ends{i} <= pos, 1);
        end

        if isempty(e)
            %!A - indel
            if (pos>=GAF.tx_start(i) && pos<=GAF.tx_stop(i)) | (pos<=GAF.tx_start(i) && pos>=GAF.tx_stop(i))
                    T.type{t} = 'Intron';
            elseif (pos<=GAF.tx_start(i) && pos<=GAF.tx_stop(i))
            		if plusstrand, T.type{t} = '5''flank';
            		else T.type{t} = '3''flank'; end
            elseif (pos>=GAF.tx_start(i) && pos>=GAF.tx_stop(i))
            		if plusstrand, T.type{t} = '3''flank';
            		else T.type{t} = '5''flank'; end
            else
            	error('Unable to map mutation to non-exonic site!');
            end
        else
            T.exon(t) = e;
            if strncmp(GAF.Gene{i},'MIR',3)         % it's a miRNA -- %%%%%NEED to fix this when actually using miRNAs in GAF
                T.type{t} = 'miRNA';
            elseif isnan(GAF.cds_start(i))        % it's a noncoding transcript
                T.type{t} = 'Non-coding_Transcript';
            elseif isempty(e_no_splice_sites)  %% it's a splice_site, promoter, or IGR(on 3' side)
                if (pos>=GAF.tx_start(i) && pos<=GAF.tx_stop(i)) | (pos<=GAF.tx_start(i) && pos>=GAF.tx_stop(i))
                    T.type{t} = 'Splice_site';
                    spliceflag = true;
                    [sd1, sd2] = deal(GAF.genomic_starts{i}(e) - pos, GAF.genomic_ends{i}(e) - pos);
                    splicedist = min(abs(sd1), abs(sd2));
                    if plusstrand; splicedist = -splicedist; end
                elseif (pos<=GAF.tx_start(i) && pos<=GAF.tx_stop(i))
            		if plusstrand, T.type{t} = '5''flank';
            		else T.type{t} = '3''flank'; end
           		elseif (pos>=GAF.tx_start(i) && pos>=GAF.tx_stop(i))
            		if plusstrand, T.type{t} = '3''flank';
            		else T.type{t} = '5''flank'; end
           		else
            		error('Unable to map mutation to non-exonic site!');
                end
            else
            
                step = 1; if ~plusstrand, step = -1; end
                g_positions = [GAF.genomic_starts{i}(e):step:GAF.genomic_ends{i}(e)];
                t_positions = [GAF.transcript_starts{i}(e):GAF.transcript_ends{i}(e)];
                cDNA_pos = t_positions(find(g_positions == pos));
                cDNA_change = change;
                if ~plusstrand, cDNA_change = rc(cDNA_change); end
                seq = F.seq{find(strcmp(T.transcript{t},F.header))};
                cDNA_old_base = seq(cDNA_pos);
                
                if cDNA_pos < GAF.cds_start(i)
                    T.type{t} = '5''-UTR';
                    seq(cDNA_pos) = cDNA_change;
                    if strfind(seq(cDNA_pos-2:cDNA_pos+2),'ATG')
                        T.type{t} = 'De_novo_Start';
                    end
                elseif cDNA_pos > GAF.cds_stop(i)
                    T.type{t} = '3''-UTR';

                else    
                        
                        mut_seq = seq;
                        mut_seq(cDNA_pos) = cDNA_change;
%
%   Note: translating the entire sequence is causing a huge slowdown to annotation....
%   We only need to translate the codon.   --ML
%
%   Also note: nt2aa is sometimes unavailable to Matlab because the toolbox license is down     --ML
%
                        translate_entire_protein = false;

                        if translate_entire_protein

                          prot_seq = nt2aa(seq(GAF.cds_start(i):GAF.cds_stop(i)),'ALTERNATIVESTARTCODONS',false);
                          mut_prot_seq = nt2aa(mut_seq(GAF.cds_start(i):GAF.cds_stop(i)),'ALTERNATIVESTARTCODONS',false);
                          
                          prot_pos = NaN;
                          for j=1:length(prot_seq)
                            if ~strcmp(prot_seq(j),mut_prot_seq(j))
                              prot_pos = j;
                              break
                            end
                          end
                          if isnan(prot_pos)
                            prot_pos = ceil((cDNA_pos - GAF.cds_start(i)+1)/3);
                          end
                        
                          [old_aa, new_aa] = deal(prot_seq(prot_pos), mut_prot_seq(prot_pos));
                          [coding_seq, mut_coding_seq] = deal(seq(GAF.cds_start(i):GAF.cds_stop(i)), mut_seq(GAF.cds_start(i):GAF.cds_stop(i)));
                          [old_codon, new_codon] = deal(coding_seq((prot_pos*3)-2:prot_pos*3), mut_coding_seq((prot_pos*3)-2:prot_pos*3));
                          [codon_start_pos, codon_end_pos] = deal(((prot_pos*3)-2)+GAF.cds_start(i)-1,(prot_pos*3)+GAF.cds_start(i)-1);

                        else

                          prot_pos = ceil((cDNA_pos - GAF.cds_start(i)+1)/3);
                          [codon_start_pos, codon_end_pos] = deal(((prot_pos*3)-2)+GAF.cds_start(i)-1,(prot_pos*3)+GAF.cds_start(i)-1);
                          [old_codon, new_codon] = deal(seq(codon_start_pos:codon_end_pos),mut_seq(codon_start_pos:codon_end_pos));
                          [old_aa, new_aa] = deal(my_nt2aa(old_codon),my_nt2aa(new_codon));
                          
                        end
                    
                end %%end of new block
                
%                if pos < R.code_start(i)
%                        if plusstrand, T.type{t} = '5''-UTR'; else T.type{t} = '3''-UTR'; end
%                        if plusstrand
%                                utr_codons = cell(3,1);
%                                %utr_starts = [pos,pos-1,pos-2]
%                                %utr_ends = [pos+2,pos+1,pos]
%                                utr_codons{1} = genome_region(chr,pos,pos+2,'hg18'); 
%                                utr_codons{2} = genome_region(chr,pos-1,pos+1,'hg18');
%                                utr_codons{3} = genome_region(chr,pos-2,pos,'hg18');
%                                for a=1:3
%                                    aa = utr_codons{a};
%                                    aa(a) = change;
%                                        if strcmp(aa,'ATG')
%                                                old_utr_codon = utr_codons{a};
%                                                utr_codon = aa;
%                                                T.type{t} = 'De_novo_Start';
%                                                cdna_pos = pos - R.tx_start(i) + 1;
%                                                T.cDNAchange{t} = sprintf('c.%d%s>%s',cdna_pos,old_base,change);
%                                                T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',cdna_pos-(a-1),cdna_pos+(3-a),old_utr_codon,utr_codon);
%                                                T.proteinchange{t} = sprintf('p.DeNovo_Start');
%                                        end
%                                end
%                        end
%                elseif pos > R.code_end(i)
%                        if plusstrand, T.type{t} = '3''-UTR'; else T.type{t} = '5''-UTR'; end
%                        if ~plusstrand
%                                utr_codons = cell(3,1);
%                                utr_codons{1} = rc(genome_region(chr,pos-2,pos,'hg18'));
%                                utr_codons{2} = rc(genome_region(chr,pos-1,pos+1,'hg18'));
%                                utr_codons{3} = rc(genome_region(chr,pos,pos+2,'hg18'));
%                                for a=1:3
%                                    aa = utr_codons{a};
%                                    aa(a) = rc(change);
%                                        if strcmp(aa,'ATG')
%                                                old_utr_codon = utr_codons{a};
%                                                utr_codon = aa;
%                                                T.type{t} = 'De_novo_Start';
%                                                cdna_pos = R.tx_end(i) - pos + 1;
%                                                T.cDNAchange{t} = sprintf('c.%d%s>%s',cdna_pos,old_base,change);
%                                                T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',cdna_pos-(a-1),cdna_pos+(3-a),old_utr_codon,utr_codon);
%                                                
%                                                T.proteinchange{t} = sprintf('p.DeNovo_Start');
%                                        end
%                                end
%                        end
%                
%                else
%                        T.type{t} = 'Coding';         % assume coding unless it's close to a splice-site
%                        spliceflag = false;
%                        cpos = pos;
%                        
%                        if e>1         % for nonfirst exons, check left splice-site
%                                jp = cpos-R.exon_starts{i}(e);
%                                if jp==-2 || jp==-1
%                                        T.type{t} = 'Splice_site';
%                                        spliceflag = true;
%                                        cpos = R.exon_starts{i}(e);         % adjust pos to within exon for subsequent annotation
%                                        splicedist = jp;
%                                end
%                        end
%                        if e<R.n_exons(i)                % for nonlast exons, check right splice-site
%                                jp = cpos-R.exon_ends{i}(e);
%                                if jp==1 || jp==2
%                                        T.type{t} = 'Splice_site';
%                                        spliceflag = true;
%                                        cpos = R.exon_ends{i}(e);         % adjust pos to within exon for subsequent annotation
%                                        splicedist = jp;
%                                end
%                        end
%                        
%                        % For coding and splice-site mutations: find position in cDNA and protein
%                        
%                        orf = [];
%                        coding_exon = 0;
%                        nframeshifts = 0;
%
%                        if plusstrand, forfrom=1; forstep=+1; forto=R.n_exons(i);
%                        else forfrom=R.n_exons(i); forstep=-1; forto=1;
%                        end
%
%                        for e=forfrom:forstep:forto
%                                st = R.exon_starts{i}(e);
%                                en = R.exon_ends{i}(e);
%                                fr = R.exon_frames{i}(e);
%                                
%        
%        %%%%%%%%%%%%%!A - what is this chunk doing???
%                                if ~isempty(orf) && fr~=-1 && mod(length(orf),3)~=fr
%                                        % recover from known programmed frameshift somewhere in the preceding exon
%                                        nframeshifts = nframeshifts+1;
%                                        pfs = mod(fr-length(orf),3);
%                                        if pfs==1, orf(end+1) = 'A';
%                                        else orf(end) = [];
%                                        end
%                                end
%                                
%        %%%%%%%%%%%%!A - builds orf and identifies orf position of mutation
%                                if st<=R.code_end(i) && en>=R.code_start(i)
%                                        coding_exon = coding_exon + 1;
%                                        if R.code_start(i)>st, st = R.code_start(i); end
%                                        if R.code_end(i)<en, en = R.code_end(i); end
%                                        if cpos >= st && cpos <= en                % mutation is within this exon
%                                                T.exon(t) = coding_exon;
%                                                if plusstrand, orf_pos = length(orf) + (cpos-st+1);
%                                                else orf_pos = length(orf) + (en-cpos+1);
%                                                end
%                                                orf_codon_start_pos = orf_pos - mod(orf_pos-1,3);
%                                                orf_codon_end_pos = orf_codon_start_pos + 2;
%                                        end
%                                        d = upper(genome_region(R.chr(i),st,en,build));
%                                        if ~plusstrand, d = rc(d); end
%                                        orf = [orf d];
%                                end
%
%                        end
%
%                        T.cDNApos(t) = orf_pos;
%                        aa_number = orf_codon_end_pos/3;
%                        T.proteinpos(t) = aa_number;
%
%                        if orf_codon_end_pos > length(orf)         % problem!        truncated ORF
%                                orf = [orf 'AA'];
%                        end
%
%                        old_codon = orf(orf_codon_start_pos:orf_codon_end_pos);
%                        if plusstrand, orf(orf_pos) = change; else orf(orf_pos) = rc(change); end
%                        new_codon = orf(orf_codon_start_pos:orf_codon_end_pos);
%                        old_aa = my_nt2aa(old_codon);
%                        new_aa = my_nt2aa(new_codon);
%
%                        if verbosity>=2
%                                T.orf{t} = orf;
%                                T.protein{t} = my_nt2aa(orf);
%                                T.nframeshifts(t) = nframeshifts;
%                        end
%
                        if spliceflag
                                T.cDNAchange{t} = sprintf('c.e%d%+d',T.exon(t),splicedist);
                                %if ~plusstrand, splicedist = -splicedist; end
                                %T.codonchange{t} = sprintf('c.e%d%+d',T.exon(t),splicedist);
                                %T.proteinchange{t} = sprintf('p.%s%d_splice',old_aa,prot_pos);
                        elseif strcmp(T.type{t}, '---') % Coding
                                if old_aa==new_aa, T.type{t} = 'Synonymous';
                                elseif new_aa=='*', T.type{t} = 'Nonsense';
                                elseif old_aa=='*', T.type{t} = 'Read-through';
                                else T.type{t} = 'Missense';
                                end

                                T.cDNAchange{t} = sprintf('c.%d%s>%s',cDNA_pos,cDNA_old_base,cDNA_change);
                                T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',codon_start_pos,codon_end_pos,old_codon,new_codon);
                                T.proteinchange{t} = sprintf('p.%s%d%s',old_aa,prot_pos,new_aa);
                        end
%                end
            end
        end
        
        
    end         % next transcript
end         

end % main function
