function T = find_mut_in_refseq(R,build,chr,st,en,change,idx,verbosity)
% find_mut_in_refseq(R, build, chromosome, start, end, change, idx, verbosity)
%
% "change" must have one of the following formats:
%        A     point mutation
%        -A    deletion
%        +CC   insertion
%        ~TT   DNP
%
% given RefSeq database in struct R (loaded using load_refseq),
% finds all Refseq transcripts overlapping the indicated mutation,
% and determines the outcome of mutation to the nucleotide specified in "change"
%
% can force annotation to a particular transcript using "idx"
%
% Mike Lawrence 2008-2010


if ~isfield(R,'gene_start'), R.gene_start = R.tx_start; end
if ~isfield(R,'gene_end'), R.gene_end = R.tx_end; end

if isnumeric(R.chr) && isnumeric(chr)
  % OK
elseif ~isnumeric(R.chr) && ~isnumeric(chr) && (strncmp(chr,'chr',3)==strncmp(R.chr{1},'chr',3))
  % OK
else   % try to reconcile incompatible chr formats
  R.chr = convert_chr(R.chr);
  chr = convert_chr(chr);
end

if ~isnumeric(st), st = str2double(st); end
if ~isnumeric(en), en = str2double(en); end

pos = st;  %variable for point mutations
is_nonSNP = false;

if ~exist('verbosity','var'), verbosity = 0; end   % flag whether to return additional info

if length(change) > 1 & (strncmp(change,'+',1) | strncmp(change,'-',1) | strncmp(change,'~',1)) %is indel or DNP
	is_nonSNP = true;
	O = struct('adjust_coordinates_for_Broad',false);  %for TCGA ovarian indel testing
	I = preprocess_nonSNP_call(R,build,chr,st,en,change,O);
	change = I.change;
	idx = I.idx;
elseif length(change)==1
  if ~exist('idx','var') || isempty(idx)
    if isnumeric(R.chr)
      chridx = find(R.chr==chr);
    else
      chridx = find(strcmp(R.chr,chr));
    end
    %if isempty(chridx), fprintf('Warning: unable to reconcile chromosome formats\n'); end
    idx = chridx(R.gene_start(chridx)<=pos & R.gene_end(chridx)>=pos);
  end
else
     error('Invalid "change" = %s',change);
end

if isempty(idx)   % IGR mutation
  idx = nan;
end

T=[];

old_base = upper(genome_region(chr,st,en,build));
change = upper(change);
non_mut_flag = strcmp(old_base,change);
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

if verbosity>=2
  T.tx_seq = repmat({'---'},length(idx),1);
  T.orf = repmat({'---'},length(idx),1);
  T.protein = repmat({'---'},length(idx),1);
  T.nframeshifts = nan(length(idx),1);
end


if is_nonSNP
	if strcmpi(I.classification(end-1:end), 'NP')  %DNP,TNP, or ONP
		I2 = find_DNP_in_refseq(R,build,I.classification,chr,I.pos1,I.pos2,change,idx,verbosity);
	elseif strcmpi(I.classification, 'DEL') | strcmpi(I.classification, 'INS')
		I2 = find_indel_in_refseq(R,build,I.classification,I.pos1,I.pos2,change,idx,verbosity);
	else
		error('Multi-nucleotide variant does not match valid types, should not happen!');
	end
	T.gene = I2.gene;
	T.strand = I2.strand;
	T.transcript = I2.transcript;
	T.type = I2.type;
	T.exon = I2.exon;
	T.cDNAchange = I2.cDNAchange;
	T.codonchange = I2.codonchange;
	T.proteinchange = I2.proteinchange;
	
	
else

for t=1:length(idx)

  i=idx(t);

  if isnan(i)
    T.type{t} = 'IGR';
    continue;
  end

  T.gene(t) = R.gene(i);
  T.strand(t) = R.strand(i);
  T.transcript(t) = R.transcript(i);

  plusstrand = strcmp(R.strand{i},'+');

  if verbosity>=2, T.tx_seq{t} = sub_get_tx_seq(i); end

  %!A - indel
  e = find(R.exon_starts{i}-2 <= pos & R.exon_ends{i}+2 >= pos, 1);

  if isempty(e)
    %!A - indel
    if pos>=R.tx_start(i) && pos<=R.tx_end(i)
      T.type{t} = 'Intron';
    else
      T.type{t} = 'Promoter';
    end

  else

    if strncmp(R.gene{i},'MIR',3)   % it's a miRNA
      T.type{t} = 'miRNA';
    elseif strncmp(R.transcript{i},'NR_',3)  % it's a noncoding transcript
      T.type{t} = 'Non-coding_Transcript';
    else
      
      if pos < R.code_start(i)
        if plusstrand, T.type{t} = '5''-UTR'; else T.type{t} = '3''-UTR'; end
        if plusstrand
          utr_codons = cell(3,1);
          %utr_starts = [pos,pos-1,pos-2]
          %utr_ends = [pos+2,pos+1,pos]
          utr_codons{1} = genome_region(chr,pos,pos+2,'hg18'); 
          utr_codons{2} = genome_region(chr,pos-1,pos+1,'hg18');
          utr_codons{3} = genome_region(chr,pos-2,pos,'hg18');
          for a=1:3
          	aa = utr_codons{a};
          	aa(a) = change;
            if strcmp(aa,'ATG')
              old_utr_codon = utr_codons{a};
              utr_codon = aa;
              T.type{t} = 'De_novo_Start';
              cdna_pos = pos - R.tx_start(i) + 1;
              T.cDNAchange{t} = sprintf('c.%d%s>%s',cdna_pos,old_base,change);
              T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',cdna_pos-(a-1),cdna_pos+(3-a),old_utr_codon,utr_codon);
              T.proteinchange{t} = sprintf('p.DeNovo_Start');
            end
          end
        end
      elseif pos > R.code_end(i)
        if plusstrand, T.type{t} = '3''-UTR'; else T.type{t} = '5''-UTR'; end
        if ~plusstrand
          utr_codons = cell(3,1);
          utr_codons{1} = rc(genome_region(chr,pos-2,pos,'hg18'));
          utr_codons{2} = rc(genome_region(chr,pos-1,pos+1,'hg18'));
          utr_codons{3} = rc(genome_region(chr,pos,pos+2,'hg18'));
          for a=1:3
          	aa = utr_codons{a};
          	aa(a) = rc(change);
            if strcmp(aa,'ATG')
              old_utr_codon = utr_codons{a};
              utr_codon = aa;
              T.type{t} = 'De_novo_Start';
              cdna_pos = R.tx_end(i) - pos + 1;
              T.cDNAchange{t} = sprintf('c.%d%s>%s',cdna_pos,old_base,change);
              T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',cdna_pos-(a-1),cdna_pos+(3-a),old_utr_codon,utr_codon);
              
              T.proteinchange{t} = sprintf('p.DeNovo_Start');
            end
          end
        end
      
      else
        T.type{t} = 'Coding';   % assume coding unless it's close to a splice-site
        spliceflag = false;
        cpos = pos;
        
        if e>1   % for nonfirst exons, check left splice-site
          jp = cpos-R.exon_starts{i}(e);
          if jp==-2 || jp==-1
            T.type{t} = 'Splice_site';
            spliceflag = true;
            cpos = R.exon_starts{i}(e);   % adjust pos to within exon for subsequent annotation
            splicedist = jp;
          end
        end
        if e<R.n_exons(i)    % for nonlast exons, check right splice-site
          jp = cpos-R.exon_ends{i}(e);
          if jp==1 || jp==2
            T.type{t} = 'Splice_site';
            spliceflag = true;
            cpos = R.exon_ends{i}(e);   % adjust pos to within exon for subsequent annotation
            splicedist = jp;
          end
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
          
		  
		  %%%%%%%%%%%%%!A - what is this chunk doing???
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
            if cpos >= st && cpos <= en    % mutation is within this exon
              T.exon(t) = coding_exon;
              if plusstrand, orf_pos = length(orf) + (cpos-st+1);
              else orf_pos = length(orf) + (en-cpos+1);
              end
              orf_codon_start_pos = orf_pos - mod(orf_pos-1,3);
              orf_codon_end_pos = orf_codon_start_pos + 2;
            end
            d = upper(genome_region(R.chr(i),st,en,build));
            if ~plusstrand, d = rc(d); end
            orf = [orf d];
          end

        end

        T.cDNApos(t) = orf_pos;
        aa_number = orf_codon_end_pos/3;
        T.proteinpos(t) = aa_number;

        if orf_codon_end_pos > length(orf)   % problem!  truncated ORF
          orf = [orf 'AA'];
        end

        old_codon = orf(orf_codon_start_pos:orf_codon_end_pos);
        if plusstrand, orf(orf_pos) = change; else orf(orf_pos) = rc(change); end
        new_codon = orf(orf_codon_start_pos:orf_codon_end_pos);
        old_aa = my_nt2aa(old_codon);
        new_aa = my_nt2aa(new_codon);

        if verbosity>=2
          T.orf{t} = orf;
          T.protein{t} = my_nt2aa(orf);
          T.nframeshifts(t) = nframeshifts;
        end

        if spliceflag
          T.cDNAchange{t} = sprintf('c.%d_splice',orf_pos);
          if ~plusstrand, splicedist = -splicedist; end
          T.codonchange{t} = sprintf('p.e%d%+d',T.exon(t),splicedist);
          T.proteinchange{t} = sprintf('p.%s%d_splice',old_aa,aa_number);
        else % Coding
          if old_aa==new_aa, T.type{t} = 'Synonymous';
          elseif new_aa=='*', T.type{t} = 'Nonsense';
          elseif old_aa=='*', T.type{t} = 'Read-through';
          else T.type{t} = 'Missense';
          end
          if ~plusstrand, T.cDNAchange{t} = sprintf('c.%d%s>%s',orf_pos,rc(old_base),rc(change));
          else T.cDNAchange{t} = sprintf('c.%d%s>%s',orf_pos,old_base,change);
          end
          T.codonchange{t} = sprintf('c.(%d-%d)%s>%s',orf_codon_start_pos,orf_codon_end_pos,old_codon,new_codon);
          T.proteinchange{t} = sprintf('p.%s%d%s',old_aa,aa_number,new_aa);
        end
      end
    end
  end
  
  if non_mut_flag, T.type{t} = 'Non-mutation'; end

end   % next transcript
end   

  % if requested, get full sequence of transcript, with exons marked by brackets and uppercase,
  % and coding sequence boundaries marked by angle-brackets.
  
  function d2 = sub_get_tx_seq(i)
    try
      d = lower(genome_region(R.chr(i),R.tx_start(i),R.tx_end(i),build));
      o = R.tx_start(i)-1;
      ins = [];  % col1 = position in d AFTER which to insert character, col2 = ascii to insert
      for e=1:R.n_exons(i)
        st = R.exon_starts{i}(e); en = R.exon_ends{i}(e);
        if en>=st
          ins = [ins;st-1 double('[');en double(']')];
          d(st-o:en-o) = upper(d(st-o:en-o));    % exons in uppercase; introns in lowercase
        end
      end
      if R.code_end(i)>=R.code_start(i)
        ins = [ins;R.code_start(i)-1 double('<');R.code_end(i) double('>')];
      end
      if size(ins,1)==0
        d2 = d;
      else
        ins = sortrows(ins);
        ins(:,1)=ins(:,1)-R.tx_start(i)+1;   % adjust coordinates so that tx_start = 1
        d2 = d(1:ins(1,1));
        for j=1:size(ins,1)
          d2 = [d2 char(ins(j,2))];
          if j<size(ins,1), d2 = [d2 d(ins(j,1)+1:ins(j+1,1))]; end
        end
        d2 = [d2 d(ins(end,1)+1:end)];
      end
    catch me
      d2 = 'error_building_transcript_sequence';
    end
  end

end % main function
