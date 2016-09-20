function run(arg1,arg2)

P=[];
P.splice_site_range = 2;

if nargin==1
  cset = arg1;
  fprintf('Assuming build hg18\n');
  build = 'hg18';
elseif nargin==2
  build = arg1;
  cset = arg2;
else 
  error('first parameter should be cset or build');
end

dir = ['/xchip/cga1/lawrence/db/' build '/effect'];
if ~exist(dir,'dir'), mkdir(dir); end

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
chr = cset(csetidx);
fprintf('\nchr%d\n',chr);
if chr<1 || chr>24, error('chr must be 1-24'); end

base = 'ACGT';

categ_list = get_effect13_categories_list;

map = nan(81,1);
map([14 32 38 40 59 65 67 33 35 39 43 47 49 15 17 23 18 24 26 36 48 52 60 62 66 70 74 76 27 63 75 79 81])=...
    [2 2 2 2 8 7 6 4 5 3 3 5 4 6 7 8 10 11 12 9 9 9 11 12 10 10 12 11 13 13 13 13 13];   % (see effect_map.xls)

if chr==1, save_struct(categ_list,[dir '/categs.txt']); end

if ~exist('loaded_flag','var')
  Rall = load_refseq(build);
  Rall.chr = convert_chr(Rall.chr);
  loaded_flag = true;
end

ln = get_chromosome_length_from_genome_file(chr,build);
R = reorder_struct(Rall,Rall.chr==chr);
R = reorder_struct(R,R.code_end>R.code_start);

E = nan(ln,4);   % Effect: nan=noncoding; 0=non-mutation; 1=change gives silent; 2=change gives nonsilent

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
      d = upper(genome_region(R.chr(i),st,en,build));
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
      if pos>ln, break; end
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
  if P.splice_site_range > 0
    for e=1:R.n_exons(i)
      if e>1
        splice_site_start = R.exon_starts{i}(e) - P.splice_site_range;
        splice_site_end = R.exon_starts{i}(e) - 1;
        E(splice_site_start:splice_site_end,:) = 2;  % non-silent
      end
      if e<R.n_exons(i)
        splice_site_start = R.exon_ends{i}(e) + 1;
        splice_site_end = R.exon_ends{i}(e) + P.splice_site_range;
        E(splice_site_start:splice_site_end,:) = 2;  % non-silent
  end,end,end

end,if i>=100, fprintf('\n'); end

% compute categ from E
E = 27*E(:,1) + 9*E(:,2) + 3*E(:,3) + E(:,4) + 1;
categ = ones(ln,1);   % 1=noncoding
idx = find(~isnan(E));
categ(idx)=map(E(idx));

% SAVE
fname = [dir '/chr' num2str(chr) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end % next chr

