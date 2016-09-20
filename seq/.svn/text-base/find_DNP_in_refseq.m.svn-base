function T = find_DNP_in_refseq(R,build,classification,chr,i_st,i_en,i_change,idx,verbosity)
% find_DNP_in_refseq(R,build,classification,i_st,i_en,change,idx,verbosity)
%
% Given Refseq database, annotates DNPs
%
%
% Alex Ramos May 2010



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

for t=1:length(idx)
  change = i_change;
  i=idx(t);

  if isnan(i)
    T.type{t} = 'IGR';
    continue;
  end
  
  T.classification{t} = classification;
  T.change{t} = change;
  T.gene(t) = R.gene(i);
  T.strand(t) = R.strand(i);
  T.transcript(t) = R.transcript(i);
  T.cDNApos_st(t) = 0;
  T.cDNApos_en(t) = 0;
  T.proteinpos_st(t) = 0;
  T.proteinpos_en(t) = 0;
  
  
  plusstrand = strcmp(R.strand{i},'+');
  e_st = find(R.exon_starts{i}-2 <= i_st & R.exon_ends{i}+2 >= i_st, 1);
  e_en = find(R.exon_starts{i}-2 <= i_en & R.exon_ends{i}+2 >= i_en, 1);
  e_start = find(R.exon_starts{i}-2 <= R.code_start(i) & R.exon_ends{i}+2 >= R.code_start(i), 1);
  e_end = find(R.exon_starts{i}-2 <= R.code_end(i) & R.exon_ends{i}+2 >= R.code_end(i), 1);
  
  e = unique([e_st e_en]);
  if length(e) > 1, e = e(1); end

  if isempty(e)
  	e = find(R.exon_starts{i}-2 > i_st & R.exon_ends{i}+2 < i_en, 1);
  	if ~isempty(e)
  	  T.type{t} = ['Multi-intron_' classification];
    elseif i_st>=R.tx_start(i) && i_en<=R.tx_end(i)
      T.type{t} = 'Intron';
    else
      T.type{t} = 'Promoter';
    end

  else %DNP falls within exon

    if strncmp(R.gene{i},'MIR',3)   % it's a miRNA
      T.type{t} = 'miRNA';
    elseif strncmp(R.transcript{i},'NR_',3)  % it's a noncoding transcript
      T.type{t} = 'Non-coding_Transcript';
    else
      cpos = [i_st i_en];
      if all([cpos < R.code_start(i)])
        if plusstrand, T.type{t} = '5''-UTR'; else T.type{t} = '3''-UTR'; end
      elseif all([cpos > R.code_end(i)])
        if plusstrand, T.type{t} = '3''-UTR'; else T.type{t} = '5''-UTR'; end
      else
        T.type{t} = 'Coding';   % assume coding unless it's close to a splice-site
        spliceflag = false;
		init_met_stop_codon_flag = false;
        
        if e == e_start | e == e_end
        	ap = cpos - R.code_start(i);
        	bp = cpos - R.code_end(i);
          	if any([ap >= 0 & ap <= 2]) | (ap(1) < 0 & ap(2) > 2) %DNP at init. Methionine
      			if plusstrand
                    iout = 'InitMet';
					T.type{t} = ['Init_Met_' classification];
                else
                	iout = 'StopCodon';
					T.type{t} = ['Stop_Codon_' classification];
                end
      			cpos([cpos < R.code_start(i)]) = R.code_start(i);  %trim
                init_met_stop_codon_flag = true;
      		elseif any([bp >= -2 & bp <= 0]) | (bp(1) < -2 & bp(2) > 0)
      			if plusstrand
                    iout = 'StopCodon';
					T.type{t} = ['Stop_Codon_' classification];
                else
                    iout = 'InitMet';
					T.type{t} = ['Init_Met_' classification];
                end
      			cpos([cpos > R.code_end(i)]) = R.code_end(i);  %trim
                init_met_stop_codon_flag = true;
        	end
        end
        
        if (i_en - i_st + 1) ~= length(change)
			error('"i_en - i_st + 1" should equal number of nucleotides in change variable');
		elseif (cpos(2)-cpos(1) + 1) ~= length(change)  % length of nucleotides in change variable differs from length of nucleotides in orf
			ii = (i_st:i_en)';
			change = change(find(cpos(1) == ii(:)):find(cpos(2) == ii(:)));  %trim change variable
		end
		
        if e>1   % for nonfirst exons, check left splice-site
          jp = cpos-R.exon_starts{i}(e);
          if (any([jp==-2 | jp==-1]) & all([jp < 0])) | (any([jp<0]) & any([jp>=0]))
            T.type{t} = ['Splice_Site_' classification];
            spliceflag = true;
            cpos = [R.exon_starts{i}(e) R.exon_starts{i}(e)];   % adjust pos to within exon for subsequent annotation
            if max(jp) == 0, splicedist = -1; else splicedist = max(jp); end
          end
          cpos([cpos < R.exon_starts{i}(e)]) = R.exon_starts{i}(e); %trim
        end
        if e<R.n_exons(i)    % for nonlast exons, check right splice-site
          jp = cpos-R.exon_ends{i}(e);
          if any([jp==2 | jp==1]) & all([jp > 0]) | (any([jp>0]) & any([jp<=0]))
            T.type{t} = ['Splice_Site_' classification];
            spliceflag = true;
            cpos = [R.exon_ends{i}(e) R.exon_ends{i}(e)];   % adjust pos to within exon for subsequent annotation
            splicedist = min(jp);
          end
          cpos([cpos > R.exon_ends{i}(e)]) = R.exon_ends{i}(e);  %trim
        end

        
        orf = [];
        coding_exon = 0;
        nframeshifts = 0;

        if plusstrand, forfrom=1; forstep=+1; forto=R.n_exons(i);
        else forfrom=R.n_exons(i); forstep=-1; forto=1;
        end

        for e=forfrom:forstep:forto
          st = R.exon_starts{i}(e);
          en = R.exon_ends{i}(e);
          fr = R.exon_frames{i}(e);
          
          if ~isempty(orf) && fr~=-1 && mod(length(orf),3)~=fr
            % recover from known programmed frameshift somewhere in the preceding exon
            nframeshifts = nframeshifts+1;
            pfs = mod(fr-length(orf),3);
            if pfs==1, orf(end+1) = 'A';
            else orf(end) = [];
            end
          end
          
		  %%%%%%%%%%%%!A - builds orf and identifies orf position of mutation
          if st<=R.code_end(i) && en>=R.code_start(i)
            coding_exon = coding_exon + 1;
            if R.code_start(i)>st, st = R.code_start(i); end
            if R.code_end(i)<en, en = R.code_end(i); end
            if all([cpos >= st & cpos <= en])  % indel is within this exon
              T.exon(t) = coding_exon;
              if plusstrand
				orf_pos = length(orf) + (cpos-st+1);
              else 
				orf_pos = length(orf) + (en-cpos+1);
				orf_pos = fliplr(orf_pos);
              end
              orf_codon_start_pos = orf_pos - mod(orf_pos-1,3);
              orf_codon_end_pos = orf_codon_start_pos + 2;
            end
            d = upper(genome_region(R.chr(i),st,en,build));
            if ~plusstrand, d = rc(d); end
            orf = [orf d];
          end

        end
        
        T.cDNApos_st(t) = orf_pos(1);
        T.cDNApos_en(t) = orf_pos(2);
        aa_number_st = orf_codon_end_pos(1)/3;
        aa_number_en = orf_codon_end_pos(2)/3;
        T.proteinpos_st(t) = aa_number_st;
        T.proteinpos_en(t) = aa_number_en;

        if orf_codon_end_pos(2) > length(orf)   % problem!  truncated ORF
          orf = [orf 'AA'];
        end



        old_codon = orf(orf_codon_start_pos(1):orf_codon_end_pos(2));
        old_aa = my_nt2aa(old_codon);
        old_bases = orf(orf_pos(1):orf_pos(2));
        
        if ~spliceflag %& ~init_met_stop_codon_flag
        	if plusstrand, orf(orf_pos(1):orf_pos(2)) = change; else, orf(orf_pos(1):orf_pos(2)) = rc(change); end
        	new_codon = orf(orf_codon_start_pos(1):orf_codon_end_pos(2));
        	new_aa = my_nt2aa(new_codon);
        end
        	
        if spliceflag
          T.cDNAchange{t} = sprintf('c.%d_splice',orf_pos(1));
          if ~plusstrand, splicedist = -splicedist; end
          T.codonchange{t} = sprintf('c.e%d%+d',T.exon(t),splicedist);
          T.proteinchange{t} = sprintf('p.%s%d_splice',old_aa,aa_number_st);
        elseif init_met_stop_codon_flag
            T.cDNAchange{t} = sprintf('c.%d_%s',orf_pos(1),iout);
            T.codonchange{t} = sprintf('c.e%d_%s',T.exon(t),iout);
            T.proteinchange{t} = sprintf('p.%s%d%s',old_aa,aa_number_st,new_aa);
        else
         	if length(orf(orf_pos(1):orf_pos(2))) ~= length(change), error('length of change different from length of ORF bases being mutated'); end
         	
         	if plusstrand, orf(orf_pos(1):orf_pos(2)) = change; else, orf(orf_pos(1):orf_pos(2)) = rc(change); end
        	new_codon = orf(orf_codon_start_pos(1):orf_codon_end_pos(2));
        	new_aa = my_nt2aa(new_codon);
         
         
          	if strcmp(old_aa,new_aa), T.type{t} = 'Synonymous';
          	elseif strcmp(new_aa(new_aa=='*'),'*'), T.type{t} = 'Nonsense';
          	elseif strcmp(old_aa(old_aa=='*'),'*'), T.type{t} = 'Read-through';
          	else T.type{t} = 'Missense';
          	end
          	
          	if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%d%s>%s',orf_pos(1),orf_pos(2),rc(old_bases),rc(change));
          	else T.cDNAchange{t} = sprintf('c.%d_%d%s>%s',orf_pos(1),orf_pos(2),old_bases,change);
          	end
          	T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',orf_codon_start_pos(1),orf_codon_end_pos(2),old_codon,new_codon);
          	
          	if aa_number_st == aa_number_en
          		T.proteinchange{t} = sprintf('p.%s%d%s',old_aa,aa_number_st,new_aa);
          	else
          		T.proteinchange{t} = sprintf('p.%s%d>%s',old_aa,aa_number_st,new_aa);
          	end

        end
      end
    end
  end
  

end   % next transcript
end
