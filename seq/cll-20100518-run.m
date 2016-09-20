
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PATIENT LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = '/xchip/cga1/lawrence/cll/analysis/20100518/freeze/CLL_paper_freeze/Individual_Set/CLL_WGS';
[a b] = system(['ls -1 ' dr '/*/wgs/coding/mut/*.maf.annotated']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = tolines(b); c = grepv('No match',c);
wgs = parse(c,[dr '/(.*)/wgs/coding/mut/.*\.maf.annotated'],{'name'});
wgs.maf = c;
wgs.wig = regexprep(wgs.maf,'.maf.annotated$','.coverage.wig.txt');
wgs.genomic_maf = regexprep(wgs.maf,'coding/','');
wgs.genomic_wig = regexprep(wgs.wig,'coding/','');


% CAPTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = '/xchip/cga1/lawrence/cll/analysis/20100518/freeze/CLL_paper_freeze/Individual_Set/CLL_Capture';
[a b] = system(['ls -1 ' dr '/*/capture/mut/*.maf.annotated']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = tolines(b); c = grepv('No match',c);
cap = parse(c,[dr '/(.*)/capture/mut/.*\.maf.annotated'],{'name'});
cap.maf = c;
cap.wig = regexprep(cap.maf,'.maf.annotated$','.coverage.wig.txt');

pat = [];
pat.name = union(cap.name,wgs.name);
pat.wgs_wig = map_across(pat.name,wgs.name,wgs.wig);
pat.wgs_maf = map_across(pat.name,wgs.name,wgs.maf);
pat.wgs_genomic_wig = map_across(pat.name,wgs.name,wgs.genomic_wig);
pat.wgs_genomic_maf = map_across(pat.name,wgs.name,wgs.genomic_maf);
pat.cap_wig = map_across(pat.name,cap.name,cap.wig);
pat.cap_maf = map_across(pat.name,cap.name,cap.maf);
pat.type = ~cellfun('isempty',pat.wgs_maf)+1;

try
  demand_file([cap.wig;cap.maf;wgs.wig;wgs.maf;wgs.genomic_maf;wgs.genomic_wig]);
  anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_all';
  save_struct(pat,[anstem '.patients.txt']);
catch me
  error('Some needed files don''t exist');
end

% make patient set of just WGS patients
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_all';
pat = load_struct([anstem '.patients.txt']);
pat = reorder_struct(pat,~cellfun('isempty',pat.wgs_wig));
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs';
save_struct(pat,[anstem '.patients.txt']);

% make patient set of just the two primary-tumor WGS individuals
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_all';
pat = load_struct([anstem '.patients.txt']);
pat = reorder_struct(pat,~cellfun('isempty',pat.wgs_wig));
pat = reorder_struct(pat,grepv('08|TT',pat.name,1));
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
save_struct(pat,[anstem '.patients.txt']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   ANALYSIS OF NONCODING MUTATIONS
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get total coverage
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
pat = load_struct([anstem '.patients.txt']);
demand_fields(pat,{'name'}); demand_file(pat.wgs_genomic_wig);
clear g; g = Genome;
for i=1:slength(pat)
  g.loadWiggle(pat.wgs_genomic_wig{i},'sum');
  fprintf('(%d/%d) Loaded %s: total coverage = %f\n',i,slength(pat),pat.name{i},g.totContents);
end
g.saveWiggle([anstem '.samples_covered.wig']);

% get WGS mutations and categories
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
pat = load_struct([anstem '.patients.txt']);
x=cell(slength(pat),1);
for i=1:slength(pat)
  x{i} = load_struct(pat.wgs_genomic_maf{i});
  x{i}.patient = repmat({pat.name{i}},slength(x{i}),1);
end
x = concat_structs(x); x.dataset = repmat({'WGS'},slength(x),1); x = add_simple_fieldnames(x);
x.chr = convert_chr(x.chr); x = make_numeric(x,{'start','end'});
x.zgr29 = get_context(x.chr,x.start,'/xchip/cga1/lawrence/db/zgr29');
x.gcsz29p = get_context(x.chr,x.start,'/xchip/cga1/lawrence/db/gcsz29p');
save_struct(x,[anstem '.genomic.maf']);

% significance analysis
regfile = '/xchip/cga1/lawrence/mm/analysis/20100318_cons/regulatory_regions_29mammal.mat';
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
mutfile = [anstem '.genomic.maf'];
covfile = [anstem '.samples_covered.wig'];
outstem = [anstem '.consreg'];
mutrate = 1.5e-6;
analyze_genomic_mutations_on_conserved_regions(regfile,mutfile,covfile,mutrate,outstem);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  mutation rate analysis using 65,536-category system
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract global coverge from cbb
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
pat = load_struct([anstem '.patients.txt']);
dr = '/xchip/cga1/lawrence/cll/analysis/20100518/freeze/CLL_paper_freeze/Individual';
tumdir = regexprep(pat.name,'(.*)',[dr '/$1/wgs/zonecov/$1-Tumor.covbb']);
normdir = regexprep(tumdir,'Tumor','Normal');
demand_file([tumdir;normdir]);
covdir = [anstem '.globcov.gcsz29p']; ensure_dir_exists(covdir);
outfile = regexprep(pat.name,'(.*)',[covdir '/$1.global_coverage.gcsz29p.txt']);
P=[]; P.covbb_file_extension = 'txt';
P.mincateg=1; P.maxcateg=65537;
get_global_coverage_stats([tumdir normdir],outfile,'gcsz29p',P);

% gather coverage stats
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
pat = load_struct([anstem '.patients.txt']);
covdir = [anstem '.globcov.gcsz29p'];
for i=1:slength(pat)
  tmp = load_struct([covdir '/' pat.name{i} '.global_coverage.gcsz29p.txt'],'%f%s%f%f%f%f%f%f%f%f%f');
  if i==1, Z=tmp; else
    Z.terr=Z.terr+tmp.terr; Z.tseqbp=Z.tseqbp+tmp.tseqbp;
    Z.nseqbp=Z.nseqbp+tmp.nseqbp; Z.tseqreads=Z.tseqreads+tmp.tseqreads;
    Z.nseqreads=Z.nseqreads+tmp.nseqreads; Z.callablebp = Z.callablebp+tmp.callablebp;
end,end

% add mutation counts
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
x = load_struct([anstem '.genomic.maf']); % (gcsz29p field added above)
Z.nmuts = histc(x.gcsz29p,1:65537);
Z.nmuts(end+1) = sum(Z.nmuts);
save_struct(Z,[anstem '.global_stats.gcsz29p.txt']);

% analyze results
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
infile = [anstem '.global_stats.gcsz29p.txt'];
outfile = [anstem '.mutrates_table.gcsz29p.txt'];
analyze_gcsz29p_stats(infile,outfile);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  WRCY hotspot analysis
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get genomic distribution of distances to nearest WRCY hotspot
dr = '/xchip/cga1/lawrence/db/';
for c=1:24, disp(c);
  load([dr '/wrcy/chr' num2str(c) '.mat'],'categ');
  load([dr '/goodbad/chr' num2str(c) '.goodbad.mat'],'f');
  minlen=min(length(categ),length(f));
  categ=categ(1:minlen); f=f(1:minlen);
  tmp = histc(categ(f==1),0:1001);
  if c==1, h=tmp; else h=h+tmp; end
end
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll-crc-mm-v2';
save([anstem '.reference.wrcy_dist.mat'],'h');

% for CLL WGS mutations, get distribution of distances to nearest WRCY hotspot
% compare to CRC, MM, PR, and MEL

anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/cll_wgs2';
cll = load_struct([anstem '.genomic.maf']);
cll.ttype = repmat({'CLL'},slength(cll),1);

anstem = '/xchip/cga1/lawrence/crc/analysis/20100524/crc9';
crc = load_struct([anstem '.genomic.maf']);
crc.ttype = repmat({'CRC'},slength(crc),1);

anstem = '/xchip/cga1/lawrence/mm/analysis/20100520/mm_wgs23';
mm = load_struct([anstem '.whole_genome_mutations.maf']);
mm.ttype = repmat({'MM'},slength(mm),1);
mm = rename_field(mm,'indiv','patient');

anstem = '/xchip/cga1/lawrence/pr/analysis/20100615/pr7';
pr = load_struct([anstem '.highconf_only.genomic.maf']);
pr.ttype = repmat({'PR'},slength(pr),1);

anstem = '/xchip/cga1/lawrence/mel/analysis/20100601/mel5wgs';
mel = load_struct([anstem '.genomic.maf']);
mel.ttype = repmat({'MEL'},slength(mel),1);
mel = rename_field(mel,'indiv','patient');

x = concat_structs_keep_all_fields({cll,crc,mm,pr,mel});
count(x.patient,1);

x.wrcy = get_context(x.chr,x.start,'/xchip/cga1/lawrence/db/wrcy');
x.goodbad = get_context(x.chr,x.start,'/xchip/cga1/lawrence/db/goodbad','.goodbad.mat');
x.goodbad(isnan(x.goodbad))=1;
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/hotspot_analysis_v3';
save_struct(x,[anstem '.genomic.maf']);
[u ui uj] = unique(x.patient);
h = nan(1002,length(u));
for i=1:length(u), disp(i)
  idx = find(uj==i & x.goodbad==1);
  h(:,i) = histc(x.wrcy(idx),0:1001);
end
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/hotspot_analysis_v3';
save([anstem '.genomic.maf.wrcy_dist.mat'],'h','u');

% combine distributions
anstem = '/xchip/cga1/lawrence/cll/analysis/20100518/hotspot_analysis_v3';
obs = load([anstem '.genomic.maf.wrcy_dist.mat']);
exp = load([anstem '.reference.wrcy_dist.mat']);
h = [exp.h obs.h]; u = ['genome';obs.u];
save_matrix(h,[anstem '.obs_exp.wrcy_dist.txt']);

% plot distributions

h2 = cumsum(bsxfun(@rdivide,h,sum(h,1)),1);
linecolor = 0.8*ones(length(u),3);
linewidth = ones(length(u),1);
idx = grep('genome',u,1); linewidth(idx) = 4; linecolor(idx,:) = [0 0 0];
idx = grep('CLL-JE',u,1); linewidth(idx) = 4; linecolor(idx,:) = [0.2 0.7 0];
idx = grep('JN_97',u,1); linewidth(idx) = 4; linecolor(idx,:) = [1 1 0];
idx = grep('MM',u,1); linecolor(idx,:) = repmat([0 0.5 1],length(idx),1);
idx = grep('CRC',u,1); linecolor(idx,:) = repmat([0.7 0 1],length(idx),1);
idx = grep('PR',u,1); linecolor(idx,:) = repmat([1 1 1],length(idx),1);
idx = grep('MEL',u,1); linecolor(idx,:) = repmat([1 0 0],length(idx),1);

maxdist = 18;
toplot = [grep('genome',u,1)];
toplot = [grep('CRC|MM',u,1);grep('genome',u,1)];
toplot = [grep('CRC|MM',u,1);grep('JN|JE|genome',u,1)];
toplot = [grep('CRC|MM',u,1);grep('JN|JE|genome',u,1);grep('PR|MEL',u,1)];
figure(1);clf;hold on;set(gca,'color',[0.8 0.8 0.8]);
for i=1:length(toplot), j=toplot(i);
  plot(0:maxdist,h2(1:maxdist+1,j),'linewidth',linewidth(j),'color',linecolor(j,:));
end
hold off, xlabel('distance (bp) to nearest AID hotspot (WRCY/RGYW)','fontsize',20);
set(gca,'fontsize',20); ylabel('fraction of mutations','fontsize',20);ylim([0 0.92]);

% MEL-0009 shows about the same level of enrichment for WRCY proximity...
%    is this effect due simply to number of mutations?

figure(2)
samplepoint = 18;
nmuts = histc(uj+1,1:length(u));
scatter(nmuts(2:end)/1000,h2(samplepoint+1,2:end),30,linecolor(2:end,:),'filled');
set(gca,'color',[0.8 0.8 0.8]);
xlabel('number of mutations (thousands)','fontsize',20);
set(gca,'fontsize',20); ylabel('frac muts within 18bp of hotspot','fontsize',20);ylim([0 0.92]);
ylim([0.78 0.93]);
