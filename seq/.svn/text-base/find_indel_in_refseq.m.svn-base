function T = find_indel_in_refseq(R,build,classification,i_st,i_en,i_change,idx,verbosity)
% find_indel_in_refseq(R,build,classification,i_st,i_en,change,idx,verbosity)
%
% Given Refseq database, annotates indels
%
%
% Alex Ramos May 2010


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
  	  T.type{t} = 'Multi-intron_deletion';
    elseif i_st>=R.tx_start(i) && i_en<=R.tx_end(i)
      T.type{t} = 'Intron';
    else
      T.type{t} = 'Promoter';
    end

  else %indel falls within exon

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
          	if any([ap >= 0 & ap <= 2]) | (ap(1) < 0 & ap(2) > 2) %indel at init. Methionine
      			if plusstrand
                    iout = 'InitMet';
                    if is_del
                        T.type{t} = 'Init_Met_Del';
                    else
                        T.type{t} = 'Init_Met_Ins';
                    end
                else
                    iout = 'StopCodon';
                    if is_del
                        T.type{t} = 'Stop_Codon_Del';
                    else
                        T.type{t} = 'Stop_Codon_Ins';
                    end
                end
      			cpos([cpos < R.code_start(i)]) = R.code_start(i);  %trim
                init_met_stop_codon_flag = true;
      		elseif any([bp >= -2 & bp <= 0]) | (bp(1) < -2 & bp(2) > 0)
      			if plusstrand
                    iout = 'StopCodon';
                    if is_del
                        T.type{t} = 'Stop_Codon_Del';
                    else
                        T.type{t} = 'Stop_Codon_Ins';
                    end
                else
                    iout = 'InitMet';
                    if is_del
                        T.type{t} = 'Init_Met_Del';
                    else
                        T.type{t} = 'Init_Met_Ins';
                    end
                end
      			cpos([cpos > R.code_end(i)]) = R.code_end(i);  %trim
                init_met_stop_codon_flag = true;
        	end
        end
		
		%if (i_en - i_st + 1) ~= length(change)
		%	error('"i_en - i_st + 1" should equal number of nucleotides in change variable');
		%elseif (cpos(2)-cpos(1) + 1) ~= length(change)  % length of nucleotides in change variable differs from length of nucleotides in orf
		%	ii = (i_st:i_en)';
		%	change = change(find(cpos(1) == ii(:)):find(cpos(2) == ii(:)));  %trim change variable
		%end
		
        if e>1   % for nonfirst exons, check left splice-site
          jp = cpos-R.exon_starts{i}(e);
          if (any([jp==-2 | jp==-1]) & all([jp < 0])) | (any([jp<0]) & any([jp>=0]))
            if is_del, T.type{t} = 'Splice_Site_Del'; else T.type{t} = 'Splice_Site_Ins'; end
            spliceflag = true;
            cpos = [R.exon_starts{i}(e) R.exon_starts{i}(e)];   % adjust pos to within exon for subsequent annotation
            if max(jp) == 0, splicedist = -1; else splicedist = max(jp); end
          end
          cpos([cpos < R.exon_starts{i}(e)]) = R.exon_starts{i}(e); %trim
        end
        if e<R.n_exons(i)    % for nonlast exons, check right splice-site
          jp = cpos-R.exon_ends{i}(e);
          if any([jp==2 | jp==1]) & all([jp > 0]) | (any([jp>0]) & any([jp<=0]))
            if is_del, T.type{t} = 'Splice_Site_Del'; else T.type{t} = 'Splice_Site_Ins'; end
            spliceflag = true;
            cpos = [R.exon_ends{i}(e) R.exon_ends{i}(e)];   % adjust pos to within exon for subsequent annotation
            splicedist = min(jp);
          end
          cpos([cpos > R.exon_ends{i}(e)]) = R.exon_ends{i}(e);  %trim
        end


        
        T.frameness{t} = '---';
        if mod(length([cpos(1):cpos(2)]),3) ~= 0 & is_del %frameshift deletion
        	is_inframe = false;
        	T.frameness{t} = 'out';
        elseif mod(length(change),3) ~= 0 & ~is_del  %frameshift insertion
        	is_inframe = false;
        	T.frameness{t} = 'out';
        else
        	is_inframe = true;
        	T.frameness{t} = 'in';
        end
        % For coding and splice-site mutations: find position in cDNA and protein
        
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

        if verbosity>=2
          T.orf{t} = orf;
          T.protein{t} = my_nt2aa(orf);
          T.nframeshifts(t) = nframeshifts;
        end


		
        old_codon = orf(orf_codon_start_pos(1):orf_codon_end_pos(2));
        old_aa = my_nt2aa(old_codon);
        if ~spliceflag
        	if is_inframe
        		even_codon_deletion = false;	
        		if is_del
        			c_del = orf(orf_pos(1):orf_pos(2));
        			%if ~plusstrand, c_del = rc(c_del); end
					new_codon = regexprep(old_codon, c_del, '');
					if strcmp(new_codon, ''), even_codon_deletion = true; end
				elseif ~is_del  %insertion
	  				if orf_codon_start_pos(1) == orf_codon_start_pos(2)  %insertion falls within a codon
						if ~plusstrand, change1 = rc(change); else change1 = change; end		
    					if orf_pos(1) - orf_codon_start_pos(1) == 1  %insertion falls within 2nd&3rd bases of codon
        					new_codon = [old_codon(1:2) change1 old_codon(end)];
        				elseif orf_pos(1) - orf_codon_start_pos(1) == 0  %insertion falls within 1st&2nd bases of codon
        					new_codon = [old_codon(1) change1 old_codon(2:end)];
        				end
        				tween_codons = false;
        			else  %insertion falls between two codons
        				tween_codons = true;
        				if plusstrand, new_codon = change; else new_codon = rc(change); end
        				%next_aa = my_nt2aa(orf(orf_codon_start_pos(1)+3:orf_codon_end_pos(2)+3));
        				next_aa = old_aa(end);
        				next_aa_pos = aa_number_st + 1;
       				end
				end
        		new_aa = my_nt2aa(new_codon);
        	end
        end
        
        
        	
        if spliceflag
          T.cDNAchange{t} = sprintf('c.%d_splice',orf_pos(1));
          if ~plusstrand, splicedist = -splicedist; end
          T.codonchange{t} = sprintf('c.e%d%+d',T.exon(t),splicedist);
          T.proteinchange{t} = sprintf('p.%s%d_splice',old_aa,aa_number_st);
        elseif init_met_stop_codon_flag
            T.cDNAchange{t} = sprintf('c.%d_%s',orf_pos(1),iout);
            T.codonchange{t} = sprintf('c.e%d_%s',T.exon(t),iout);
            T.proteinchange{t} = sprintf('p.%s%d_%s',old_aa,aa_number_st,iout);
        elseif ~is_inframe
 			if strcmpi(classification, 'ins')
 				T.type{t} = 'Frame_Shift_Ins';
 				if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%dins%s',orf_pos(1),orf_pos(2),rc(change));
          		else T.cDNAchange{t} = sprintf('c.%d_%dins%s',orf_pos(1),orf_pos(2),change);
          		end
 			elseif strcmpi(classification, 'del')
 				T.type{t} = 'Frame_Shift_Del';
 				%if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),rc(change));
          		%else T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),change);
          		%end
				if  ~plusstrand cc = rc(change);, else cc = change; end
				if orf_pos(1) == orf_pos(2)
					T.cDNAchange{t} = sprintf('c.%ddel%s',orf_pos(1),cc);
				else
					T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),cc);
				end
				
        	end
        	T.codonchange{t} = sprintf('c.(%d-%d)%sfs',orf_codon_start_pos(1),orf_codon_end_pos(2),old_codon);
			T.proteinchange{t} = sprintf('p.%s%dfs',old_aa(1),aa_number_st);
        elseif is_inframe % Coding
			if ~is_del
				T.type{t} = 'In_Frame_Ins';
				if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%dins%s' ,orf_pos(1),orf_pos(2),rc(change));
				else T.cDNAchange{t} = sprintf('c.%d_%dins%s',orf_pos(1),orf_pos(2),change);
				end
				
				if tween_codons
					T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',orf_codon_start_pos(1),orf_codon_end_pos(2),old_codon,[old_codon(1:3) new_codon old_codon(4:end)]);
					T.proteinchange{t} = sprintf('p.%s%d_%s%dins%s',old_aa(1),aa_number_st,next_aa(1),next_aa_pos,new_aa);
				else
					T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',orf_codon_start_pos(1),orf_codon_end_pos(2),old_codon,new_codon);
					T.proteinchange{t} = sprintf('p.%s%d>%s',old_aa,aa_number_st,new_aa);
				end
			elseif is_del
				T.type{t} = 'In_Frame_Del';
				%if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),rc(c_del));
				%else T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),c_del);
				%end
				if orf_pos(1) == orf_pos(2)
					T.cDNAchange{t} = sprintf('c.%ddel%s',orf_pos(1),c_del);
				else
					T.cDNAchange{t} = sprintf('c.%d_%ddel%s',orf_pos(1),orf_pos(2),c_del);
				end
				
				
				if even_codon_deletion
					T.codonchange{t} = sprintf('c.(%d-%d)%sdel',orf_codon_start_pos(1),orf_codon_end_pos(2),old_codon);
					if length(old_aa) > 1
						T.proteinchange{t} = sprintf('p.%s%d_%s%ddel',old_aa(1),aa_number_st,old_aa(end),aa_number_en);
					else
						T.proteinchange{t} = sprintf('p.%s%ddel',old_aa,aa_number_st);
					end
				elseif length(old_aa) > 1
					T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',orf_codon_start_pos(1),orf_codon_end_pos(2),old_codon,new_codon);
					T.proteinchange{t} = sprintf('p.%s%d_%s%d>%s',old_aa(1),aa_number_st,old_aa(end),aa_number_en,new_aa);
				else
					T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',orf_codon_start_pos(1),orf_codon_end_pos(2),old_codon,new_codon);
					T.proteinchange{t} = sprintf('p.%s%d>%s',old_aa,aa_number_st,new_aa);
				end		
			end
        end
      end
    end
  end
  

end   % next transcript
end
