function [RG miR]=make_hg_file(hg_output_fname,refgene_fname,reflink_fname, ...
                          refseqstatus_fname,wg_fname,add_miR_to_RG, ...
                               miR_output_fname)

% from http://genome.ucsc.edu/goldenPath/gbdDescriptions.html#GenePredictions
%
% Gene Predictions and RefSeq Genes
% The following definition is used for gene prediction tables. 
% In alternative-splicing situations, each transcript has a row in this table. 
% table genePred
% "A gene prediction."
%     (
%     string  name;               "Name of gene"
%     string  chrom;              "Chromosome name"
%     char[1] strand;             "+ or - for strand"
%     uint    txStart;            "Transcription start position"
%     uint    txEnd;              "Transcription end position"
%     uint    cdsStart;           "Coding region start"
%     uint    cdsEnd;             "Coding region end"
%     uint    exonCount;          "Number of exons"
%     uint[exonCount] exonStarts; "Exon start positions"
%     uint[exonCount] exonEnds;   "Exon end positions"
%     )


%rg=read_table(refgene_fname,'%s%s%s%d%d%d%d%d%s%s',char(9),0);
% It seems they added an id as first column
%rg=read_table(refgene_fname,'%s%s%s%s%d%d%d%d%d%s%s%d%s%s%s%s',char(9),0);
% It seems they moved back to the prev format
if (0)
  extra_col=0;
  rg=read_table(refgene_fname,'%s%s%s%d%d%d%d%d%s%s',char(9),0);
else
  extra_col=1;
  rg=read_table(refgene_fname,'%s%s%s%s%d%d%d%d%d%s%s%d%s%s%s%s',char(9),0);
end
rg_ln=line_count(refgene_fname);
if size(rg.dat{1},1)~=rg_ln
  error('did not read all lines in refgene file');
end

% Adding miRNA data if wg table present (download from http://hgdownload.cse.ucsc.edu/goldenPath/hg...
if exist('wg_fname','var') & ~isempty(wg_fname)
  wg=read_dlm_file(wg_fname);
  wg_ln=line_count(wg_fname);
  if size(wg,2)~=wg_ln
    error('did not read all lines in wg file')
  end
  count = 0;
  for i = 1:wg_ln
    if strmatch('mirna',lower(wg{i}(10)))
      count = count+1;
      miR(count).bin = str2num(wg{i}{1});
      miR(count).chr = wg{i}{2};
      miR(count).start = str2num(wg{i}{3});
      miR(count).end = str2num(wg{i}{4});
      miR(count).name = wg{i}{5};
      miR(count).score = str2num(wg{i}{6});
      if wg{i}{7} == '+'
        miR(count).strand = 1;
      elseif wg{i}{7} == '-'
        miR(count).strand = 0;
      else
        miR(count).strand = wg{i}{7};
        warning(['strand for miR ' num2str(count) ' unexpected: ' miR(count).strand])
      end
      miR(count).thick_start=str2num(wg{i}{8});
      miR(count).thick_end=str2num(wg{i}{9});
    end
  end
else
  miR = [];
end


% from http://genome.ucsc.edu/goldenPath/gbdDescriptions.html#GenePredictions
%
% RefSeq Link
% First used Dec. 2000 
% table refLink
% "Link between RefSeq mRNAs and HUGO, LocusLink etc."
%     (
%     string name;        "Name displayed in UI"
%     string product;     "Name of protein product"
%     string mrnaAcc;     "mRNA accession"
%     string protAcc;     "Protein accession"
%     uint geneId;        "Pointer to geneName table"
%     uint prodId;        "Pointer to prodName table"
%     uint locusLinkId;   "Locus Link ID"
%     uint omimId;        "OMIM ID"
%     )

rl=read_table(reflink_fname,'%s%s%s%s%d%d%d%d',char(9),0);
rl_ln=line_count(reflink_fname);
if size(rl.dat{1},1)~=rl_ln
  error('did not read all lines');
end

% from http://genome.ucsc.edu/goldenPath/gbdDescriptions.html#GenePredictions
%
%RefSeq Status
%First used Dec. 2001 
%table refSeqStatus
%"Links RefSeq mRNA accessions with status"
%    (
%    string mrnaAcc;     "RefSeq mRNA accession"
%    string status;      "RefSeq status (Reviewed, Provisional, Predicted)"
%    )

if exist('refseqstatus_fname','var')
  rs=read_table(refseqstatus_fname,'%s%s%s',char(9),0);
  rs_ln=line_count(reflink_fname);
  if size(rs.dat{1},1)~=rs_ln
    error('did not read all lines');
  end
end

[nms,rli]=sort(rl.dat{3});
[Ms,rsi,ms2]=match_string_sets_hash(rs.dat{1},nms);
if range(ms2-(1:max(ms2))') ~=0
  error('should match');
end

[Mg,rgi,mg2]=match_string_sets_hash(rg.dat{1+extra_col},nms);

num_genes = length(mg2);

for j=1:num_genes
  i=mg2(j);
  if mod(j,1000)==0
    disp(j);
  end
  RG(j).refseq=nms{i};
  RG(j).gene=rl.dat{2}{rli(i)};
  RG(j).symb=rl.dat{1}{rli(i)};
  RG(j).locus_id=rl.dat{7}(rli(i));
  RG(j).chr=rg.dat{2+extra_col}{rgi(j)};
  RG(j).strand=convert_enum(rg.dat{3+extra_col}{rgi(j)},{'+',1;'-',0});
  RG(j).start=rg.dat{4+extra_col}(rgi(j));
  RG(j).end=rg.dat{5+extra_col}(rgi(j));
  RG(j).cds_start=rg.dat{6+extra_col}(rgi(j));
  RG(j).cds_end=rg.dat{7+extra_col}(rgi(j));
  RG(j).status=rs.dat{2}{rsi(i)};
end

if add_miR_to_RG % toggle to combine miRNA data to RG, despite not all
                 % terms being present. Locus_id defined to -1*bin.
  for j = 1: length(miR)
    RG(num_genes+j).refseq = miR(j).name;
    RG(num_genes+j).gene = miR(j).name;
    RG(num_genes+j).symb = miR(j).name;
    RG(num_genes+j).locus_id = -1*miR(j).bin;
    RG(num_genes+j).chr = miR(j).chr;
    RG(num_genes+j).strand = miR(j).strand;
    RG(num_genes+j).start = miR(j).start;
    RG(num_genes+j).end = miR(j).end;
    RG(num_genes+j).cds_start = miR(j).thick_start;
    RG(num_genes+j).cds_end = miR(j).thick_end;
    RG(num_genes+j).status = 'miRNA';
  end
end

f=fopen(hg_output_fname,'w');
fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','refseq','gene','locuslink','chro','strand','start','end');
for i=1:length(RG);
  if mod(i,1000)==0
    disp(i);
  end
  fprintf(f,'%s\t%s\t%d\t%s\t',RG(i).refseq,RG(i).gene,RG(i).locus_id,RG(i).chr);
  if RG(i).strand==1
    fprintf(f,'+\t');
  else
    fprintf(f,'-\t');
  end
  fprintf(f,'%d\t%d\n',RG(i).start,RG(i).end);
end
fclose(f);

if exist('wg_fname','var') & ~isempty(wg_fname) & ~add_miR_to_RG
  f=fopen(miR_output_fname,'w');
  fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\n','gene','bin','chro','strand','start','end');
  for i=1:length(miR);
    if mod(i,1000)==0
      disp(i);
    end
    fprintf(f,'%s\t%d\t%s\t',miR(i).name,miR(i).bin,miR(i).chr);
    if miR(i).strand==1
      fprintf(f,'+\t');
    else
      fprintf(f,'-\t');
    end
    fprintf(f,'%d\t%d\n',miR(i).start,miR(i).end);
  end
  fclose(f);
end
