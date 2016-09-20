function run(arg1,arg2)

trackname = 'effect29b';
write_text_file = false;

P=[];
%P.splice_site_range = 2;
P.splice_site_range_intron = 2;    % number of bases on each side of intron to designate "splice_site"
P.splice_site_range_exon = 2;      % number of bases on each side of exon to designate "splice_site"

if nargin==1
  cset = arg1;
elseif nargin==2
  build = arg1;
  cset = arg2;
end

if ~exist('build','var')
  fprintf('Assuming build hg19\n');
  build = 'hg19';
end

dir = ['/xchip/cga1/lawrence/db/' build '/' trackname];
if ~exist(dir,'dir'), mkdir(dir); end

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
chr = cset(csetidx);
fprintf('\nchr%d\n',chr);
if chr<1 || chr>24, error('chr must be 1-24'); end

base = 'ACGT';

categ_list = get_effect29_categories_list;

if chr==1, save_struct(categ_list,[dir '/categs.txt']); end

if ~exist('loaded_flag','var')
  Rall = load_refseq(build);
  Rall.chr = convert_chr(Rall.chr);
  loaded_flag = true;
end

ln = get_chromosome_length_from_genome_file(chr,build);
R = reorder_struct(Rall,Rall.chr==chr);
R = reorder_struct(R,R.code_end>R.code_start);

E = zeros(ln,3);   % Effect: 1=change gives silent; 2=change gives missense 3=change gives nonsense
                   %         4=splice-site                    (note: nonstop is treated as nonsense)

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

  orfi = listmap(orf,base);
  
  % TRY ALL MUTATIONS
  for c=1:3:3*floor(length(orf)/3)  % for each codon
    old_codon = orf(c:c+2);
    old_aa = my_nt2aa(old_codon);
    for j=1:3  % for each position in the codon
      orf_pos = c+j-1;
      pos = genome_pos(orf_pos);
      if pos>ln, break; end
      if plusstrand, genome_refi = orfi(orf_pos); else genome_refi = 5-orfi(orf_pos); end
      genome_altis = setdiff(1:4,genome_refi);
      for k=1:3     % for each possible substitution
        genome_alti = genome_altis(k);
        if plusstrand, orf_alti = genome_alti; else orf_alti = 5-genome_alti; end
        new_codon = old_codon; new_codon(j) = base(orf_alti); new_aa = my_nt2aa(new_codon);
        if old_aa~=new_aa
          if new_aa=='*'
            E(pos,k) = max(E(pos,k),3);           % nonsense
          elseif old_aa=='*'
            E(pos,k) = max(E(pos,k),3);           % nonstop (treated as nonsense)
          else
            E(pos,k) = max(E(pos,k),2);           % missense
          end
        else
          E(pos,k) = max(E(pos,k),1);             % synonymous
  end,end,end,end

  % mark splice-sites
  if P.splice_site_range_intron > 0
    for e=1:R.n_exons(i)
      if e>1
        splice_site_start = R.exon_starts{i}(e) - P.splice_site_range_intron;
        splice_site_end = R.exon_starts{i}(e) - 1;
        E(splice_site_start:splice_site_end,:) = 4;  % splice
      end
      if e<R.n_exons(i)
        splice_site_start = R.exon_ends{i}(e) + 1;
        splice_site_end = R.exon_ends{i}(e) + P.splice_site_range_intron;
        E(splice_site_start:splice_site_end,:) = 4;  % splice
  end,end,end

  if P.splice_site_range_exon > 0
    for e=1:R.n_exons(i)
      splice_site_start = R.exon_starts{i}(e) + 0;
      splice_site_end = R.exon_starts{i}(e) + (P.splice_site_range_exon-1);
      E(splice_site_start:splice_site_end,:) = 4;  % splice

      splice_site_start = R.exon_ends{i}(e) - (P.splice_site_range_exon-1);
      splice_site_end = R.exon_ends{i}(e) + 0;
      E(splice_site_start:splice_site_end,:) = 4;  % splice
  end,end

end,if i>=100, fprintf('\n'); end

% compute categ from E
categ = 9*(E(:,1)-1) + 3*(E(:,2)-1) + 1*(E(:,3)-1) + 1;
categ(E(:,1)==4) = 28;  % splice-site
categ(E(:,1)==0) = 29;  % noncoding

% SAVE
fname = [dir '/chr' num2str(chr) '.mat'];
save(fname,'categ');

if write_text_file
  f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
  fprintf(f,'%d\n',categ);
  fclose(f);
end

end % next chr

