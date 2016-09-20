function C = load_cosmic_database(P)
% Mike Lawrence 2008-2012

if ~exist('P','var'), P=[]; end

if ischar(P)
  if P(1)=='/'
    fprintf('Assuming "%s" is the COSMIC file to load.\n',P);
    tmp = P;
    P=[];
    P.cosmic_file = tmp;
  else     
    fprintf('Assuming "%s" is the build.\n',P);
    tmp = P;
    P=[];
    P.build = tmp;
  end
end

if ~(isfield(P,'build') && ~isempty(P.build)) && ~(isfield(P,'build_cosmic') && isempty(P.build_cosmic))
  fprintf('Assuming hg18\n');
  P.build = 'hg18';
end
P = impose_default_value(P,'build_cosmic',P.build);
P = impose_default_value(P,'cosmic_file',[]);
P = impose_default_value(P,'length_filter',10);

% load and prep COSMIC database

if strcmp(P.build_cosmic,'hg18_old')
  P = impose_default_value(P,'cosmic_file','/cga/tcga-gsc/home/lawrence/xchip_tcga_gbm_analysis_lawrence/db/CosmicMutantExport_v43_260809.tsv');

  % STEPS FOR DOWNLOADING DATABASE
  %
  % (1) wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicMutantExport_v43_260809.tsv
  % (2) edit header line to be tab-delimited (was no longer necessary as of v42)
  %     in v43, the first two names on the header line were space-delimited: had to manually edit to tab
  % (3) remove extra blank line between header line and first data line
  %     -- still necessary as of v43

  C = load_struct(P.cosmic_file);
  if isfield(C,'NCBI36genomestop')
    C = rename_field(C,'NCBI36genomestop','end');
    C = rename_field(C,'NCBI36genomestart','start');
    C = rename_field(C,'Chromosome','chr');
  elseif isfield(C,'NCBI36genomeposition')
    C = parse_in(C,'NCBI36genomeposition','([^:]*):(\d*)-(\d*)',{'chr','start','end'});
  else
    error('Can''t find genomic coordinates fields in %s!',P.cosmic_file);
  end
else   % all other versions

  %% wget ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicMutantExport_v62_291112.tsv.gz

  P = impose_default_value(P,'cosmic_file','/cga/tcga-gsc/home/lawrence/xchip_tcga_gbm_analysis_lawrence/db/CosmicMutantExport_v62_291112.tsv');
  C = load_struct(P.cosmic_file);
  if isfield(C,'NCBI36genomeposition')
    C = rename_field(C,'NCBI36genomeposition','MutationNCBI36genomeposition');
  end
  if isfield(C,'GRCh37genomeposition')
    C = rename_field(C,'GRCh37genomeposition','MutationGRCh37genomeposition');
  end
  if strcmp(P.build_cosmic,'hg18') || strcmp(P.build_cosmic,'hg18_v2')
    demand_field(C,'MutationNCBI36genomeposition');

    if isfield(C,'MutationGRCh37genomeposition')
      % hack for missing hg18 coordinates of XPO1 mutation  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      idx = grep('2\:61719472-61719472',C.MutationGRCh37genomeposition,1);
      C.MutationNCBI36genomeposition(idx) = repmat({'chr2:61572976-61572976'},length(idx),1);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    C = parse_in(C,'MutationNCBI36genomeposition','([^:]*):(\d*)-(\d*)',{'chr','start','end'});

  elseif strcmp(P.build_cosmic,'hg19')
    demand_field(C,'MutationGRCh37genomeposition')
    C = parse_in(C,'MutationGRCh37genomeposition','([^:]*):(\d*)-(\d*)',{'chr','start','end'});
  else
    error('Unknown P.build_cosmic %s',P.build_cosmic);
  end
end

%keyboard
fprintf('Loaded COSMIC database from %s\n', P.cosmic_file);
C = rename_field(C,'Genename','gene');
C.chr = convert_chr(C.chr);
C = make_numeric(C,{'start','end'});

% at this point, there is the opportunity to rescue mutations
%   that have protein coordinates but not DNA coordinates
C = parse_in(C,'MutationAA','^p.(\D)(\d+)(\D)',{'aa_old','aa_start','aa_new'},2);

goodidx = find(~isnan(C.start) & ~isnan(C.aa_start));
badidx = find(isnan(C.start) & ~isnan(C.aa_start));
[u ui uj] = unique(C.gene(badidx));
for i=1:length(u)
  gidx = goodidx(strcmp(u{i},C.gene(goodidx)));
  bidx = badidx(uj==i);
  z = listmap(C.aa_start(bidx),C.aa_start(gidx));
  m = find(~isnan(z));
  if ~isempty(m) % matched some!
    C.chr(bidx(m)) = C.chr(z(m));
    C.start(bidx(m)) = C.start(z(m));
    C.end(bidx(m)) = C.end(z(m));
  end
end
% (next step would be to map to transcripts)



C = reorder_struct(C,~isnan(C.chr)&~isnan(C.start)&~isnan(C.end));
len = C.end-C.start+1;
C = reorder_struct(C,len<=P.length_filter);

C.site = cell(length(C.chr),1);
for i=1:length(C.chr)
  c = C.chr(i);
  s = C.start(i);
  e = C.end(i);
  if ~isempty(c) && isnumeric(s) && isnumeric(e)
    C.site{i} = sprintf('%s|chr%d:%d-%d', C.gene{i},c,s,e);
  else
    C.site{i} = '';
  end
end

C.site = regexprep(C.site,'chr23','chrX');
C.site = regexprep(C.site,'chr24','chrY');

tmp=C;
C=[];
C.mut=tmp;
[C.gene.name tmp C.mut.gene] = unique(C.mut.gene);
C.gene.n_muts = zeros(slength(C.gene),1);
for i=1:slength(C.gene)
  C.gene.n_muts(i) = sum(C.mut.gene==i);
end
