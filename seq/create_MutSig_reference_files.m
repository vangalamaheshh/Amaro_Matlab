function create_MutSig_reference_files(build_dir, P)
% create_MutSig_reference_files(build_dir, P)
%
% INPUT PARAMETER
%
% build_dir          directory where input files will be read from, and output files will be written.
%
% INPUT FILES
%
% build_info.txt,    filename of tab-separated headered textfile with the following format:
%
%         num     chromosome     length          male    female  use     alias
%         1       chr1           247249719       2       2       D       1
%         2       chr2           242951149       2       2       D       2
%         3       chr3           199501827       2       2       D       3
%         ...
%         21      chr21          46944323        2       2       D       21
%         22      chr22          49691432        2       2       D       22
%         23      chrX           154913754       1       2       D       X
%         24      chrY           57772954        1       0       D       Y
%         25      chrM           16571           999     999     X       M,MT,chrMT
%         26      chr1_random    1663265         2       2       X       1_random
%         27      chr2_random    185571          2       2       X       2_random
%         28      chr3_random    749256          2       2       X       3_random
%         ...
%
%     NOTE: This table is read by the Java class ReferenceInfo, which looks for it as:
%           build_dir = "/path/build"   -->  build_info_file = "/path/build/build_info.txt"
%     NOTE: "alias" column is optional.
%     NOTE: These names will be used to search for the reference sequence flatfiles in the
%           directory specified by the parameter flatfiles_dir
%
% chr*.txt,       reference sequence chr1.txt, chr2.txt, etc.
%
%     NOTE: These files should contain no whitespace and no characters other than acgtACGTN.
%           They should be same length (in bytes) as the lengths (in basepairs) in the build_info_file.
%
% refGene.txt,    filename of tab-separated unheadered textfile with the following columns in this order:
%         id                  1189
%         transcript          NM_173528
%         chr                 chr15
%         strand              +
%         tx_start            79213698
%         tx_end              79228571
%         code_start          79213725
%         code_end            79227929
%         n_exons             7
%         exon_starts         79213698,79214665,79215911,79217446,79223057,79227253,79227733,
%         exon_ends           79213794,79214755,79216144,79217531,79223216,79227327,79228571,
%         tmp                 0
%         gene                C15orf26
%         tmp                 cmpl
%         tmp                 cmpl
%         exon_frames         0,0,0,2,0,0,2,
%         version             1
% 
%      NOTE: This is the format of the UCSC "refGene.txt" file made available for each UCSC build,
%            from e.g. wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/refGene.txt.gz
%      NOTE: "version" column is optional.
%
%
% OUTPUT FILES
%
% R.mat             Matlab MAT file version of refgene_table
%
% target_list.txt   list of refGene exons, tab-delimited, unheadered:
%
%         OR4F5   1       69089   70010   4.273319e-01    922     1
%         OR4F16  1       367657  368597  4.601488e-01    941     1
%         OR4F29  1       367657  368597  4.601488e-01    941     1
%         OR4F3   1       367657  368597  4.601488e-01    941     1
%
%      NOTE: columns = gene, chr#, start, end, frac_GC, length, mem
%      NOTE: "mem" (membership) column is no longer used; OK to set to all 1
%      NOTE: start and end are both expanded outward by 2bp to include splice-sites.
%
% genelist.txt      list of genes in target_list.txt
%
% all.fwb           trinculeotide context and mutation effect at each genomic basepair    
% all.fwi           (in FWB format, with FWI index file)
%
% categs.txt        key to the values in all.fwb
%
%
% Mike Lawrence 2012-07-11

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'splice_site_bp_to_add',2);

%%%%%%%%   INPUTS
%%%%%%%%

if ~exist(build_dir','dir'), error('directory not found: %s',build_dir); end

build = regexprep(build_dir,'^.*/([^/]*)$','$1');
build_info_file = fullfile(build_dir,[build '_info.txt']);
refgene_table = fullfile(build_dir,'refGene.txt');

demand_file(build_info_file);
demand_file(refgene_table);

% READ BUILD INFO FILE

fprintf('Loading build info file\n');
ReferenceInfoObj.init(build_dir);
maxchr = ReferenceInfoObj.getMaxNum(build);
C = [];
cidx = 1;
for chr=1:maxchr
  try
    use = ReferenceInfoObj.getUse(chr,build);
    if strcmpi(use,'D')
      C.num(cidx,1) = chr;
      C.chromosome(cidx,1) = ReferenceInfoObj.getChrom(chr,build);
      C.length(cidx,1) = ReferenceInfoObj.getLength(chr,build);
      cidx=cidx+1;
    elseif ~stcmpi(use,'X')
      error('"use" column in build_info_file should contain only X or D.');
    end
  catch me
    % chromosome not defined
  end
end
nchr = cidx-1;

% VERIFY EXISTENCE OF REFERENCE SEQUENCE FLATFILES

fprintf('Verifying existence of reference sequence flatfiles\n');
C.flatfile = cell(nchr,1);
notfound = [];
for ci=1:nchr, chr=C.num(ci);
  found = false;
  name = [build_dir '/chr' num2str(chr) '.txt'];
  if exist(name,'file')
    found = true;
  else
    name = [build_dir '/chr' C.chromosome{chr} '.txt'];
    if exist(name,'file')
      found = true;
%    else
%      if isfield(B,'alias')
%        aliases = split(C.aliases{chr},',');
%        if ~isempty(aliases)
%          for i=1:length(aliases)
%            name = [build_dir '/chr' aliases{i} '.txt'];
%            if exist(name,'file')
%              found = true;
%              break
%      end,end,end,end
    end
  end
  if found
    C.flatfile{chr} = name;
  else
    notfound = [notfound chr];
end,end

if ~isempty(notfound)
  disp(notfound);
  error('Could not find flatfile(s) for the above chromosome(s)');
end

% READ REFGENE TABLE

fprintf('Loading refGene table\n');
colnames = {...
    'id','transcript','chr','strand','tx_start','tx_end','code_start','code_end',...
    'n_exons','exon_starts','exon_ends','col12','gene','col14','col15','exon_frames','version'};
R = load_struct_noheader(refgene_table,colnames);

fprintf('Converting refGene table\n');
R = make_numeric(R,{'n_exons','tx_start','code_start','tx_end','code_end'});
z = zeros(slength(R),1);
R.tx_len = z; R.code_len = z;
for i=1:slength(R)
  if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
  R.tx_start(i) = R.tx_start(i)+1;
  R.code_start(i) = R.code_start(i)+1;
  R.exon_starts{i} = str2double(split(R.exon_starts{i}(1:end-1),','))+1;
  R.exon_ends{i} = str2double(split(R.exon_ends{i}(1:end-1),','));
  R.exon_frames{i} = str2double(split(R.exon_frames{i}(1:end-1),','));
  for e=1:R.n_exons(i)
    st = R.exon_starts{i}(e);
    en = R.exon_ends{i}(e);
    R.tx_len(i) = R.tx_len(i) + (en-st+1);
    if en<R.code_start(i) || st>R.code_end(i), continue; end
    if st<R.code_start(i), st = R.code_start(i); end
    if en>R.code_end(i), en = R.code_end(i); end
    R.code_len(i) = R.code_len(i) + (en-st+1);
  end
end,fprintf('\n');
R.n_codons = R.code_len / 3;


%%%%%%%%   OUTPUTS
%%%%%%%%

% MATLAB VERSION OF REFGENE TABLE

outfile = [build_dir '/R.mat'];
if exist(outfile,'file')
  fprintf('Already exists: %s\n',outfile);
else
  fprintf('Writing %s\n',outfile);
  save(outfile,'R');
end

% TARGET LIST

outfile = [build_dir '/target_list.txt'];
if exist(outfile,'file')
  fprintf('Already exists: %s\n',outfile);
else
  fprintf('Generating target list\n');

  R.chr = convert_chr(R.chr,build);
  R = reorder_struct(R,~isnan(R.chr));   % remove genes on non-used chromosomes
  R = reorder_struct(R,R.code_len>0);    % remove RNA-coding genes
  
  E = cell(slength(R),1);
  for i=1:slength(R)
    st = max(R.code_start(i),R.exon_starts{i}-P.splice_site_bp_to_add);
    en = min(R.code_end(i),R.exon_ends{i}+P.splice_site_bp_to_add);
    idx = find(en<st);
    st(idx)=[]; en(idx)=[];
    if ~isempty(st)
      E{i}.gene = repmat(R.gene(i),length(st),1);
      E{i}.chr = repmat({num2str(R.chr(i))},length(st),1);  % needs to be string for condense_regions
      E{i}.start = st;
      E{i}.end = en;
    end
  end
  E = concat_structs(E);
  E = condense_regions(E);

  fprintf('Annotating GC content\n');
  E.chr = str2double(E.chr);
  E.gc = get_gc_content(E.chr,E.start,E.end,build);
  
  E.len = E.end-E.start+1;
  
  E.mem = ones(slength(E),1);

  E = sort_struct(E,{'chr','start','end'});

  fprintf('Writing %s\n',outfile);
  save_struct_noheader(E,outfile);
end

% GENE LIST

g = unique(E.gene);
outfile = [build_dir '/genelist.txt'];
if exist(outfile,'file')
  fprintf('Already exists: %s\n',outfile);
else
  fprintf('Writing %s\n',outfile);
  save_lines(g,outfile);
end

% CONTEXT65EFFECT files

% categs.txt

c65 = generate_categ_context65_names();
e13 = get_effect13_categories_list();

X=[];
X.num = (1:65*13)';
X.name = cell(slength(X),1);
x=1;
for c=1:65, for e=1:13
  X.name{x} = [c65.name{c} ':' e13.name{e}];
  x=x+1;
end,end

outfile = [build_dir '/categs.txt'];
if exist(outfile,'file')
  fprintf('Already exists: %s\n',outfile);
else
  fprintf('Writing %s\n',outfile);
  save_struct(X,outfile);
end

% all.fwi

Y=[];
Y.chr = C.num;
Y.start = ones(slength(Y),1);
Y.end = C.length;

outfile = [build_dir '/all.fwi'];
if exist(outfile,'file')
  fprintf('Already exists: %s\n',outfile);
else
  fprintf('Writing %s\n',outfile);
  save_struct_noheader(Y,outfile);
end

% all.fwb

outfile = [build_dir '/all.fwb'];
if exist(outfile,'file')
  fprintf('Already exists: %s\n',outfile);
else
  fprintf('Writing %s\n',outfile);

  Rall = R;
  base = 'ACGT';
  map = nan(81,1);   % (see /xchip/cga1/lawrence/db/hg18/effect/effect_map.xls)
  map([14 32 38 40 59 65 67 33 35 39 43 47 49 15 17 23 18 24 26 36 48 52 60 62 66 70 74 76 27 63 75 79 81])=...
      [2 2 2 2 8 7 6 4 5 3 3 5 4 6 7 8 10 11 12 9 9 9 11 12 10 10 12 11 13 13 13 13 13];

  outwidth = 16;     % in order to encompass 1-845
  endianness = 'b';  % crucial because Matlab is little-endian by default, and file will be read by Java
  outtype = ['uint' num2str(outwidth)];
  fout = fopen(outfile,'w');

  for ci=1:nchr
    chr = C.num(ci); len = C.length(ci);
    fprintf('chr%d: ',chr);

    % get context65

    dna = upper(genome_region(chr,1,len,build));
    if length(dna)<len
      dna = [dna repmat('N',len-length(dna),1)];
    end

    val = zeros(1,len);
    val(dna=='A') = 1;
    val(dna=='C') = 2;
    val(dna=='G') = 3;
    val(dna=='T') = 4;

    context65 = 16*(val(2:end-1)-1) + 4*(val(1:end-2)-1) + val(3:end);
    context65(val(2:end-1)==0 | val(1:end-2)==0 | val(3:end)==0) = 65;
    context65 = [65 context65 65]';

    % get effect13

    R = reorder_struct(Rall,Rall.chr==chr);
    R = reorder_struct(R,R.code_end>R.code_start);

    E = nan(len,4);   % Effect: nan=noncoding; 0=non-mutation; 1=change gives silent; 2=change gives nonsilent

    for i=1:slength(R)
      if ~mod(i,100), fprintf('%d/%d ',i,slength(R)); end

      % BUILD ORF
      orf = []; genome_pos = []; nframeshifts = 0;
      plusstrand = strcmp(R.strand{i},'+');
      if plusstrand, forfrom=1; forstep=+1; forto=R.n_exons(i);
      else forfrom=R.n_exons(i); forstep=-1; forto=1;
      end
      for e=forfrom:forstep:forto
        st = R.exon_starts{i}(e); en = R.exon_ends{i}(e); fr = R.exon_frames{i}(e);
        % recover from known programmed frameshift somewhere in the preceding exon
        if ~isempty(orf) && fr~=-1 && mod(length(orf),3)~=fr
          nframeshifts = nframeshifts+1; pfs = mod(fr-length(orf),3);
          if pfs==1
            orf(end-1:end) = [];
            genome_pos(end-1:end) = [];
          else
            orf(end) = [];
            genome_pos(end) = [];
          end
        end
        % look up sequence and add to orf
        if st<=R.code_end(i) && en>=R.code_start(i)
          if R.code_start(i)>st, st = R.code_start(i); end
          if R.code_end(i)<en, en = R.code_end(i); end
          if st>=1 && en<=len
            d = dna(st:en);
          elseif st>=1 && en>len
            d = [dna(st:len) repmat('N',en-len,1)];
          elseif st>len
            d = repmat('N',en-st+1,1);
          end
          if plusstrand, p=st:en; else p=en:-1:st; d = rc(d); end
          orf = [orf d]; genome_pos = [genome_pos p];
        end
      end

      % TRY ALL MUTATIONS
      for c=1:3:3*floor(length(orf)/3)  % for each codon
        old_codon = orf(c:c+2);
        old_aa = my_nt2aa(old_codon);
        for j=1:3  % for each position in the codon
          pos = genome_pos(c+j-1);
          if pos>len, break; end
          for k=1:4  % for each possible substitution
            if plusstrand, genomebase=k; else genomebase=5-k; end
            if old_codon(j)==base(k)
              E(pos,genomebase) = 0;       % non-mutation
            else
              new_codon = old_codon; new_codon(j) = base(k); new_aa = my_nt2aa(new_codon);
              if old_aa~=new_aa
                E(pos,genomebase) = 2;     % non-silent
              else
                if E(pos,genomebase)~=2    % (if no other non-silent change found)
                  E(pos,genomebase) = 1;   % silent
      end,end,end,end,end,end
    
      % mark splice-sites as "any change is nonsilent"
      if P.splice_site_bp_to_add > 0
        for e=1:R.n_exons(i)
          if e>1
            splice_site_start = R.exon_starts{i}(e) - P.splice_site_bp_to_add;
            splice_site_end = R.exon_starts{i}(e) - 1;
            E(splice_site_start:splice_site_end,:) = 2;  % non-silent
          end
          if e<R.n_exons(i)
            splice_site_start = R.exon_ends{i}(e) + 1;
            splice_site_end = R.exon_ends{i}(e) + P.splice_site_bp_to_add;
            E(splice_site_start:splice_site_end,:) = 2;  % non-silent
      end,end,end

    end,fprintf('\n');

    % compute effect from E
    E = 27*E(:,1) + 9*E(:,2) + 3*E(:,3) + E(:,4) + 1;
    effect = ones(len,1);   % 1=noncoding
    idx = find(~isnan(E));
    idx(idx>len) = [];
    effect(idx)=map(E(idx));

    % compute context65effect
    c65e = (context65-1)*13 + effect;

    % write to file
    fwrite(fout,c65e,outtype,endianness);
  end    % next chromosome

  fclose(fout);

end

fprintf('Done!\n');










