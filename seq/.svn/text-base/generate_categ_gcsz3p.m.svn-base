function generate_categ_gcsz3p(cset)

% creates gcsz3p category files
%
% non-N, non-end = 1-65536
%
% G  0 32768             good            bad
% C  0 16384             nonconserved    conserved
% S  0 4096 8192 12288   neither_transcribed    (+)transcribed    (-)transcribed   both_transcribed
% Z  0 1024 2048 3072    IGR     intron     UTR    exon
% P  0-1023              (pentamer_sequence)
%    +1
%
% N/end = 65537

outdir = '/xchip/cga1/lawrence/db/gcsz3p';

t = nan(65537,24);
for cidx=1:length(cset)
  chr=cset(cidx)

  cname = ['chr' num2str(chr)]; fprintf('%s\n',cname);
  fprintf('Loading\n');
  g = load(['/xchip/cga1/lawrence/db/goodbad/' cname '.goodbad.mat']); g=g.f; % T=good F=bad
  c = load(['/xchip/cga1/lawrence/db/regulatory/' cname '.mat']); c=c.r;      % F=noncons T=cons
  s = load(['/xchip/cga1/lawrence/db/strand/' cname '.mat']); s=s.strand;     % 0=nei 1=+ 2=- 3=both
  z = load(['/xchip/cga1/lawrence/db/zone/' cname '.mat']); z=z.zone;         % 0=IGR 1=int 2=UTR 3=ex
  p = load(['/xchip/cga1/lawrence/db/pentamer/' cname '.mat']); p=p.categ;    % 1-1024 / 1025=N
  
  fprintf('Computing\n');
  maxlen = max([length(g),length(c),length(s),length(z),length(p)]);
  if maxlen>length(g), g(end+1:maxlen)=0; end
  if maxlen>length(c), c(end+1:maxlen)=0; end
  if maxlen>length(s), s(end+1:maxlen)=0; end
  if maxlen>length(z), z(end+1:maxlen)=0; end
  if maxlen>length(p), p(end+1:maxlen)=1025; end

  v = 32768*double(~g) + 16384*double(c) + 4096*s + 1024*z + p;
  v(p==1025) = 65537; 
  t(:,chr)=histc(v,1:65537);

  fprintf('Saving mat\n'); outstem = [outdir '/' cname]; save([outstem '.mat'],'v');
  fprintf('Saving txt\n'); out = fopen([outstem '.txt'],'wt'); fprintf(out,'%d\n',v); fclose(out);

  if chr==1
    % create categs.txt
    
    g = load_struct('/xchip/cga1/lawrence/db/goodbad/categs.txt');
    c = load_struct('/xchip/cga1/lawrence/db/regulatory29/categs.txt');
    s = load_struct('/xchip/cga1/lawrence/db/strand/categs.txt');
    z = load_struct('/xchip/cga1/lawrence/db/zone/categs.txt');
    p = load_struct('/xchip/cga1/lawrence/db/pentamer/categs.txt');
    
    x = p.name(1:1024);
    n={};for i=1:slength(z), n{end+1} = regexprep(x,'(.*)',[z.name{i} ':$1']); end; x = cat(1,n{:});
    n={};for i=1:slength(s), n{end+1} = regexprep(x,'(.*)',[s.name{i} ':$1']); end; x = cat(1,n{:});
    n={};for i=1:slength(c), n{end+1} = regexprep(x,'(.*)',[c.name{i} ':$1']); end; x = cat(1,n{:});
    n={};for i=slength(g):-1:1, n{end+1} = regexprep(x,'(.*)',[g.name{i} ':$1']); end; x = cat(1,n{:});
    
    categ = []; categ.num = (1:65537)'; categ.name = [x;'any N'];
    save_struct(categ,[outdir '/categs.txt']);
  end

end   % next chr
