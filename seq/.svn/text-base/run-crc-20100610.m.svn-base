% (1) get patient list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PATIENT LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = '/xchip/cga1/lawrence/crc/analysis/20100610/freeze';
[a b] = system(['ls -1 ' dr '/*/coding-wgs/mut/calls/*.coverage.wig.txt']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = tolines(b); c = grepv('No match',c);
wgs = parse(c,[dr '/(.*)/coding-wgs/mut/calls/.*\.coverage.wig.txt'],{'name'});
wgs.wig = c;
wgs.maf = regexprep(wgs.wig,'calls/(.*).coverage.wig.txt$','annotated/$1.maf.annotated');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = 'null';
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
pat.cap_wig = map_across(pat.name,cap.name,cap.wig);
pat.cap_maf = map_across(pat.name,cap.name,cap.maf);
pat.genomic_wig = regexprep(pat.wgs_wig,'/coding-wgs','/wgs');
pat.genomic_maf = regexprep(pat.wgs_maf,'/coding-wgs','/wgs');

anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
save_struct(pat,[anstem '.patients.txt']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MUTATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% point mutations
pm = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/mut_coding/CRC-all.maf.annotated');
map.fr = {'5EKFOAO4','5TA9VAT9','C0245_T','C0258_T','C0294_T','C1','C3','MQE9OAS7','TV28IANZ'};
map.to = {'CRC-0002','CRC-0003','CRC-0007','CRC-0008','CRC-0009','CRC-0004','CRC-0005','CRC-0010','CRC-0006'};
pm.patient = map_across(pm.Tumor_Sample_Barcode,map.fr,map.to);

% indels
in = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/indels/somatic_coding_filtered_indels.maf.annotated');
in.patient = regexprep(in.Tumor_Sample_Barcode,'-Tumor$','');

% combine
x = concat_structs_keep_all_fields({pm,in});
x.dataset = repmat({'WGS'},slength(x),1);
x = add_simple_fieldnames(x);
x = reorder_struct(x,grepv('UTR|Intron|Promoter|miRNA|IGR',x.type,1));
x.context = get_context(x.chr,x.start,'/xchip/tcga_scratch/lawrence/db/context65');
x = collapse_adjacent_mutations(x);
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
save_struct(x,[anstem '.maf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  COVERAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract coverage from wigs
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
pat = load_struct([anstem '.patients.txt']);
targlist = ['/xchip/tcga_scratch/lawrence/capture/' ...
           'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP_noheader.txt'];
categfile = '/xchip/tcga_scratch/lawrence/db/context65/weRefseq_terr_only.wig';
mincateg = 1; maxcateg = 65;
outdir = [anstem '.covfiles']; ensure_dir_exists(outdir);
outfiles = cell(slength(pat),1);
wigfiles = cell(slength(pat),1);
for i=1:slength(pat)
  outfiles{i} = [outdir '/' pat.name{i} '.somatic_coverage.txt'];
  wigfiles{i} = {};
  if ~isempty(pat.wgs_wig{i}), wigfiles{i}{end+1} = pat.wgs_wig{i}; end
  if ~isempty(pat.cap_wig{i}), wigfiles{i}{end+1} = pat.cap_wig{i}; end
end
pat.covfile = outfiles;
save_struct(pat,[anstem '.patients.txt']);
extract_from_wig(pat.name,targlist,outfiles,categfile,mincateg,maxcateg,wigfiles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load extracted coverage and save as Matlab file
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
pat = load_struct([anstem '.patients.txt']);
C = []; C.sample = pat; C.ns = slength(C.sample);
C.file.targ = ['/xchip/tcga_scratch/lawrence/capture/' ...
           'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP_noheader.txt'];
C.file.categdir ='/xchip/tcga_scratch/lawrence/db/context65'; C.ncat = 65;
C = load_coverage(C);
save([anstem '.coverage.C.mat'],'C','-v7.3');

% automated mutation category discovery
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
load([anstem '.coverage.C.mat'],'C');
N = squeeze(sum(sum(C.cov,1),2));
n = get_context_from_maf([anstem '.maf'],65);
Nn = collapse_Nn_65_to_32([N n]);
P=[];
P.max_k = 5;
P.mutcategs_report_filename = [anstem '.mutcateg_discovery.txt'];
Ks = find_mut_categs(Nn,P);
% choose what categories to use
K = Ks{4};
K.name = {'CpG_transition';'other_C_transition';'A_transition';'transversion'};
save_struct(K,[anstem '.mutcategs.txt']);

% collapse coverage
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
load([anstem '.coverage.C.mat'],'C');
%K = load_struct('/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/categs_CpG_CG_AT_Indel.txt');
K = load_struct([anstem '.mutcategs.txt']);
C1 = collapse_coverage_categories(C,K);
C1 = compute_fraction_coverage(C1);
save([anstem '.coverage.C1.mat'],'C1','-v7.3');

% DRAW COVERAGE PLOT
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
load([anstem '.coverage.C1.mat'],'C1');
P=[]; P.sample_name_xadj = 0;
P.sample_name_fontsize = 8;P.colorbar_fontsize = 7;P.GC_plot_fontsize = 7;
figure(1);clf;snameidx = capture_covplot(C1,P);
save([anstem '.covplot_samporder.mat'],'snameidx');
print('-dpng','-r120',[anstem '.covplot.png']);
print('-dpng','-r360',[anstem '.covplot_hires.png']);

% DRAW BARGRAPHS
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
I = [];
I.maffile = [anstem '.maf'];
I.covfile = [anstem '.coverage.C1.mat'];
I.patfile = [anstem '.patients.txt'];
I.ordfile = [anstem '.covplot_samporder.mat'];
D = load_mutation_data_for_bargraph(I);
P=[]; P.plot_ypos = [0.70 0.40 0.10]; P.plot_width = 0.91; P.plot_left = 0.07;P.legend = 'best';
P.y_axis_fontsize = 7; P.ylabel_fontsize = 12; P.sample_name_fontsize = 8;
P.off_scale_fontsize = 5; P.text_xadj = -0.1; P.barwidth = 0.7;
figure(2);clf;display_mutation_bargraph2(D,P);
print('-dpng','-r120',[anstem '.bargraphs.png']);
print('-dpng','-r360',[anstem '.bargraphs_hires.png']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   SIGNIFICANCE ANALYSIS
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add "indel" and "null" categories
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
M = load_struct([anstem '.maf']);
%K = load_struct([anstem '.mutcategs.txt']);
K1 = load_struct([anstem '.mutcategs.txt']);
K2 = load_struct('/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/categs_CpG_CG_AT_Indel.txt');
K = concat_structs_keep_common_fields({K1,reorder_struct(K2,4)});  % add indel category
M.categ = assign_mut_categs(M,K);
idx = grepi('Frameshift|Nonsense|Read-through|Nonstop|Non-stop|Splice',M.type,1);
M = make_numeric(M,'categ');
M.categ(idx) = 6;    % add null category
save_struct(M,[anstem '.with_null_category.maf']);
K = concat_structs({K,reorder_struct(K,1)});
K.autoname{end} = 'null'; K.name{end} = 'null'; K.type{end} = 'indel'; K.from{end} = 'N'; K.change{end} = 'null';
save_struct(K,[anstem '.mutcategs.with_null_category.txt']);

% load final data set
P = [];
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
P.mutfile = [anstem '.with_null_category.maf'];
P.covfile = [anstem '.coverage.C1.mat'];
%P.catfile = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/categs_CpG_CG_AT_Indel_Null.txt';
P.catfile = [anstem '.mutcategs.with_null_category.txt'];
P.patlist = [anstem '.patients.txt'];
P.genelist = '/xchip/tcga_scratch/lawrence/capture/weRefseq_genelist.txt';
M = load_all_mutation_data2(P);



M = analyze_mutation_data(M,anstem,P);

% (repeat)
outstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.run1a';
M = analyze_mutation_data(M,outstem,P);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   ANALYZE WHOLE-GENOME MUTATIONS IN WINDOWS
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = load_struct(['/xchip/cga1/lawrence/crc/analysis/20100610/freeze/mut_genomic/CRC-' ...
                 'all.genomic_mutations.maf.annotated']);
x = add_simple_fieldnames(x); x.chr = convert_chr(x.chr); x = make_numeric(x,'start');
window_size = 5000;
x.q = 1e6*x.chr+round(x.start/window_size);
h = histc(x.q,1:max(x.q));
[ct q] = sort(h,'descend');








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   CIRCOS PLOTS
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% update location of dRanger results
f = load_struct('/xchip/cga1/lawrence/crc/analysis/20100524/circos/circos_files.txt');
f.dRanger_results = regexprep(f.individual,'(.*)',['/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/' ...
                    '$1.dRanger_results.somatic.txt']);
save_struct(f,'/xchip/cga1/lawrence/crc/analysis/20100610/circos/circos_files.txt');

% make CIRCOS plots
dr = '/xchip/cga1/lawrence/crc/analysis/20100610/circos';
ensure_dir_exists(dr);
f = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/circos/circos_files.txt');
demand_file([f.segfile;f.dRanger_results]);
cd(dr);
libdir = '/xchip/tcga/gdac_prod/applications/genepattern/taskLib/CreateCircosPlot.19.2433/';
cex = [libdir 'circos'];
perl_lib = ['/home/radon00/lawrence/CancerGenomeAnalysis/trunk/analysis_pipeline/genepattern/common_tools/circos/' ...
            'circos-0.52/lib'];
for i=1:slength(f), disp(i);
  iname = f.individual{i};
  drf = f.dRanger_results{i};
  seg = f.segfile{i};
  fh_CreateCircosPlot(libdir,cex,perl_lib,iname,drf,seg);
%  disp('please manually re-run the command'); keyboard
% (works on nobelium: don't need to manually re-run)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   FIGURES OF COMPLEX REARRANGEMENTS
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drfile = '/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic3_rearrs_with_fixed_annotations.mat';
D = load(drfile); D=D.x;
D = reorder_struct(D,D.score>=3);
f = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/crc9_segfile_locations.txt');

% CRC-0004

dr4 = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.segseq/CRC-0004/matfiles';
ssmat4 = [dr4 '/CRC-0004_6006F_nonmedfc_tumor_normal_aligned_paired_reads_qual20.mat'];
segmat4 = [dr4 '/CRC-0004_6006F_seg_W_400_initFP_1000_pmerge_auto.mat'];
S4 = load(ssmat4);
seg4 = load(segmat4);
D4 = reorder_struct(D,grep('CRC-0004',D.individual,1));
P=[];
P.midstrip = 0.1;
P.xon = [0 20]*1e6;
P.compression_factor = 1;
P.no_labeling = true;
%P.segfile = f.seg_file{grep('CRC-0004',f.Sample_Id,1)};
P.line_anchors = 'segseq_dots';
draw_multichr_plot(S4,seg4,D4,[8 20],P)
set(gca,'position',[0.05 0.05 0.9 0.9]);


% CRC-0006

dr6 = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.segseq/CRC-0006/matfiles';
ssmat6 = [dr6 '/CRC-0006-cluster4.mat'];
segmat6 = [dr6 '/CRC-0006-cluster4_seg_W_400_initFP_1000_pmerge_auto.mat'];
S6 = load(ssmat6);
seg6 = load(segmat6);
D6 = reorder_struct(D,grep('CRC-0006',D.individual,1));
build = 'hg18'
R = load_refseq(build);
P=[];
P.compression_factor = 2;
%P.segfile = f.seg_file{grep('CRC-0006',f.Sample_Id,1)};
% P.inset_chr = 5;  P.inset_left = 122420000;  P.inset_right = 123137000; P.refseq = R;
P.no_labeling = true;
P.xon = [0 30]*1e6;
P.line_anchors = 'segseq_dots';
draw_multichr_plot(S6,seg6,D6,[5 11],P)




% CRC-0003
% (dont have good SegSeq results for this one)
D3 = reorder_struct(D,grep('CRC-0003',D.individual,1));
P=[];
P.midstrip = 0.1;
P.compression_factor = 2;
P.xon = [0 150 0]*1e6;
%P.no_labeling = true;
P.segfile = f.seg_file{grep('CRC-0003',f.Sample_Id,1)};
P.line_anchors = 'segfile';
draw_multichr_plot([],[],D3,[1 3 5],P)


% to plot SegSeq dots using Derek's function
addpath ~/CancerGenomeAnalysis/trunk/segseq/
plotRatios(S.RATIOS,seg.SEG)
set(gca,'xticklabel',[],'yticklabel',[]);xlabel('');ylabel('')


% overlap between dR and CN breakpoints

drfile = '/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic3_rearrs_with_fixed_annotations.mat';
D = load(drfile); D=D.x;
D = reorder_struct(D,D.score>=3);
f = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/crc9_segfile_locations.txt');
D.cnb1 = false(slength(D),1); D.cnb2 = D.cnb1;
margin = 100000;
for p=2:10
  name = sprintf('CRC-%04d',p);
  didx = grep(name,D.individual,1);
  s = load_cnseg_file(f.seg_file{grep(name,f.Sample_Id,1)});
  s = reorder_struct(s,s.nprobes>=8);
  for j=1:length(didx), i=didx(j);
    D.cnb1(i) = any(s.chr==D.chr1(i) & abs(s.start-D.pos1(i))<=margin | abs(s.end-D.pos1(i))<=margin);
    D.cnb2(i) = any(s.chr==D.chr2(i) & abs(s.start-D.pos2(i))<=margin | abs(s.end-D.pos2(i))<=margin);
  end
  fprintf('%s\t%d rearrs\t%d segs\n',name,length(didx),slength(s));
  D.cnb = nansub({'neither','one','both'},D.cnb1+D.cnb2+1);
  xcount(D.class(didx),D.cnb(didx));
  xcount(D.class([didx;didx]),[D.cnb1(didx);D.cnb2(didx)]);
end




% stacking plot
xcount(D.class,D.individual)
[uc uci ucj] = unique(D.class);
[ui uii uij] = unique(D.individual);
h = hist2d_fast(ucj,uij);
bar(h','stacked');
legend(uc,'interpreter','none','location','northwest')
fs=12;
ylabel('number of rearrangements','fontsize',15);
set(gca,'xticklabel',ui,'fontsize',fs)

% by location
D.z1 = regexprep(D.site1,'.*(IGR|Intron|Exon|UTR|Promoter).*','$1');
D.z2 = regexprep(D.site2,'.*(IGR|Intron|Exon|UTR|Promoter).*','$1');
D.z1 = regexprep(D.z1,'(Exon|Intron|UTR|Promoter)','Gene');
D.z2 = regexprep(D.z2,'(Exon|Intron|UTR|Promoter)','Gene');
for i=1:slength(D)
  q = sort([D.z1(i);D.z2(i)]);
  D.zz{i,1} = [q{1} '-' q{2}];
end
xcount(D.T_lenhomology,D.zz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   CRC-0009   TCF7l2-VTI1a fusion figure
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drfile = '/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic3_rearrs_with_fixed_annotations.mat';
D = load(drfile); D=D.x;
D = reorder_struct(D,D.score>=3);
look(D,490);
pos1 = 114220870;
pos2 = 114760547;

% dRanger reads
w = load('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/CRC-0009/wgs/ra/CRC-0009-Tumor.all.weird.pairs.mat'); w=w.X;
idx = find(w(:,2)==10 & abs(w(:,4)-pos1)<500 & w(:,6)==10 & abs(w(:,8)-pos2)<500);   % 13 tumreads
w = w(idx,2:end);
worig = w;
% order randomly
rand('twister',54);
w = worig(randperm(size(worig,1)),:);
% draw figure
clf;hold on;
win1left = pos1-300; win1right = pos1+50;
win2left = pos2-50; win2right = pos2+300;
leftgap = 100; midgap = 200; rightgap = 100;
win1start = leftgap; win1width = win1right-win1left;
win2start = leftgap + win1width + midgap; win2width = win2right - win2left;
grey = ones(1,3)*0.8;
rectangle('position',[win1start+win1width -1 midgap size(w,1)*2+5],'facecolor',grey);
for i=1:size(w,1)
  line(win1start+[w(i,3) w(i,4)]-win1left,i*2*[1 1],'linewidth',6,'color',[0.2 0.2 0.8]);
  line([win1start+w(i,4)-win1left win2start+w(i,7)-win2left],i*2*[1 1],'linewidth',2,'color',[0.3 0.3 0.3]);
  line(win2start+[w(i,7) w(i,8)]-win2left,i*2*[1 1],'linewidth',6,'color',[0.8 0.2 0.2]);
end
xt = [win1start win1start+win1width win2start win2start+win2width];
xtl = num2cellstr([win1left win1right win2left win2right]);
set(gca,'xtick',xt,'xticklabel',xtl,'tickdir','out','ytick',[]);
ylim([0 size(w,1)*2+2]);
hold off  



% BreakPointer reads
% fusion sequence
fseq = load_textfile('/xchip/cga1/lawrence/crc/analysis/20100610/VTI1a_TCF7L2_rearragement_fusion.txt');
d1 = upper(genome_region(10,pos1-200,pos1+200));
d2 = upper(genome_region(10,pos2-200,pos2+200));
[fseq(1:101);d1(100:200)],[fseq(102:202);d2(200:300)]
firstwin1pos = pos1-101;
firstwin2base = 102;
firstwin2pos = pos2-1;
% extracted fusion-spanning reads from BreakPointer file:
%     grep 14_17_10 CRC-0009-Tumor.aligned.sam
rorig = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/VTI1a_TCF7L2_rearragement_reads.txt');
rorig = make_numeric(rorig,'pos');
rorig.readlen = cellfun('length',rorig.seq);
[u ui uj] = unique(rorig.pos); rorig = reorder_struct(rorig,ui);
% order randomly
rand('twister',54);
r = reorder_struct(rorig,randperm(slength(rorig)));
% draw figure
clf;hold on
for i=1:slength(r)
  line([r.pos(i) firstwin2base],i*2*[1 1],'linewidth',6,'color',[0.2 0.2 0.8]);
  line([firstwin2base r.pos(i)+r.readlen(i)-1],i*2*[1 1],'linewidth',6,'color',[0.8 0.2 0.2]);
end
ylim([0 slength(r)*2+2]);
lf = min(r.pos); rt = max(r.pos+r.readlen);
xt = [lf firstwin2base rt];
xtl = {num2str(firstwin1pos+lf-1),[num2str(pos1) ' / ' num2str(pos2)],num2str(firstwin2pos+rt-firstwin2base)}
set(gca,'xtick',xt,'xticklabel',xtl,'tickdir','out','ytick',[]);
hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   ANALYSIS OF WHOLE-GENOME MUTATIONS
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combine WGS mutations into single file
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
pat = load_struct([anstem '.patients.txt']);
x=cell(slength(pat),1);
for i=1:slength(pat)
  x{i} = load_struct(pat.genomic_maf{i});
  x{i}.patient = repmat({pat.name{i}},slength(x{i}),1);
end
x = concat_structs(x); x.dataset = repmat({'WGS'},slength(x),1); x = add_simple_fieldnames(x);
x.chr = convert_chr(x.chr); x = make_numeric(x,{'start','end'});
save_struct(x,[anstem '.genomic.maf']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   RP REGION ANALYSIS
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get total coverage
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
pat = load_struct([anstem '.patients.txt']);
demand_fields(pat,{'name'});
clear g; g = Genome;
for i=1:slength(pat)
  g.loadWiggle(pat.genomic_wig{i},'sum');
  fprintf('(%d/%d) Loaded %s: total coverage = %f\n',i,slength(pat),pat.name{i},g.totContents);
end
g.saveWiggle([anstem '.samples_covered.wig']);

% significance analysis
regfile = '/xchip/cga1/lawrence/mm/analysis/20100318_cons/regulatory_regions_29mammal.mat';
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
mutfile = [anstem '.genomic.maf'];
covfile = [anstem '.samples_covered.wig'];
outstem = [anstem '.consreg'];
mutrate = 5e-6;
analyze_genomic_mutations_on_conserved_regions(regfile,mutfile,covfile,mutrate,outstem);



% analysis on beta-catenin binding regions
R = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/chipseq/bcatenin_regions.txt'); % Bottomly et al.
R.chr = convert_chr(R.chr);
R = make_numeric(R,{'start','end'});
save('/xchip/cga1/lawrence/crc/analysis/20100610/chipseq/bcatenin_regions.mat','R');
% significance analysis
regfile = '/xchip/cga1/lawrence/crc/analysis/20100610/chipseq/bcatenin_regions.mat';
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
mutfile = [anstem '.genomic.maf'];
covfile = [];  % (will assume full coverage)
outstem = [anstem '.bcatenin'];
mutrate = 5.6e-6;
analyze_genomic_mutations_on_conserved_regions(regfile,mutfile,covfile,mutrate,outstem);


% analysis on tcf4 binding regions
R = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/chipseq/tcf4_regions.txt'); % Hatzis et al.
R.chr = convert_chr(R.Chromosome);
R = make_numeric(R,{'TCF4peakcenter'});
R.start = R.TCF4peakcenter-300;
R.end = R.TCF4peakcenter+299;
R.gene = R.Closestgenename;
R.zone = mapacross(R.Classificationofbindingsite,...
  {'3''-proximal','5''-proximal','TSS-3''','enhancer','intragenic','unclassified'},...
  {'UTR','UTR','UTR','UTR','IGR','IGR'});
R.ig = nan(slength(R),1);
save('/xchip/cga1/lawrence/crc/analysis/20100610/chipseq/tcf4_regions.mat','R');
R2 = keep_fields(R,{'chr','start','end'});
save_struct(R2,'/xchip/cga1/lawrence/crc/analysis/20100610/chipseq/tcf4_regions_600bp.txt');
% significance analysis
regfile = '/xchip/cga1/lawrence/crc/analysis/20100610/chipseq/tcf4_regions.mat';
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
mutfile = [anstem '.genomic.maf'];
covfile = [];  % (will assume full coverage)
outstem = [anstem '.tcf4'];
mutrate = 5.6e-6;
analyze_genomic_mutations_on_conserved_regions(regfile,mutfile,covfile,mutrate,outstem);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  whole-genome: analyze mutation rates, Ka/Ks ratio, CpG islands
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract global coverge from cbb
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
pat = load_struct([anstem '.patients.txt']);
dr = '/xchip/cga1/lawrence/crc/analysis/20100610/freeze';
pat.tumcbb = regexprep(pat.name,'(.*)',[dr '/$1/wgs/zonecov/$1-Tumor.covbb']);
pat.normcbb = regexprep(pat.tumcbb,'Tumor','Normal');
demand_file([pat.tumcbb;pat.normcbb]);
save_struct(pat,[anstem '.patients.txt']);
covdir = [anstem '.globcov.' categdir]; ensure_dir_exists(covdir);
outfile = regexprep(pat.name,'(.*)',[covdir '/$1.global_coverage.' categdir '.txt']);
P=[]; P.covbb_file_extension = 'txt';
get_global_coverage_stats([pat.tumcbb pat.normcbb],outfile,categdir,P);

% annotate genomic mutations by category
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
x = load_struct([anstem '.genomic.maf']);
x = setfield(x,categdir,get_context(x.chr,x.start,['/xchip/cga1/lawrence/db/' categdir]));
save_struct(x,[anstem '.genomic.maf']);

% load per-sample coverage stats
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
pat = load_struct([anstem '.patients.txt']);
covdir = [anstem '.globcov.' categdir];
N = cell(slength(pat),1);
for i=1:slength(pat), fprintf('%d/%d ',i,slength(pat));
  N{i} = load_struct_specify_string_cols([covdir '/' pat.name{i} '.global_coverage.' categdir '.txt'],2);
end,fprintf('\n');
Ntot = 0; for i=1:length(N), Ntot=Ntot+N{i}.callablebp(1:end-1); end
save([anstem '.globcov.' categdir '.N.mat'],'N','Ntot');

% analyze Ka/Ks in pooled data
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
load([anstem '.globcov.' categdir '.N.mat'],'N','Ntot');
P=[];
x = load_struct([anstem '.genomic.maf']);
ntot = get_context_from_mafstruct(x,[],categdir);
Nn = [Ntot ntot];
outstem = [anstem '.globcov.' categdir '.analysis'];
[R H] = analyze_KaKs_data2(Nn,categdir,outstem,P);

% check "effect"
Z = get_categs(categdir);
x = make_numeric(x,categdir);
x.categname = nansub(Z.name,getfield(x,categdir));
x2 = reorder_struct(x,grepv('noncoding',x.categname,1));
x2.silent = regexprep(x2.categname,'.*change to (.*) is silent.*','$1');
x2.silent = regexprep(x2.silent,'.*any change is silent.*','ACGT');
x2.silent = regexprep(x2.silent,'.*any change is nonsilent.*','none');
for i=1:slength(x2), x2.sbs(i,1) = ismember(x2.tum_allele1{i},x2.silent{i})|ismember(x2.tum_allele2{i},x2.silent{i});end
x3=reorder_struct(x2,grepv('bad',x2.categname,1))
x3=reorder_struct(x3,~cellfun('isempty',x3.categname));
xcount(x3.type,x3.sbs);
idx = find(strcmp(x3.type,'Synonymous')&x3.sbs==0);


% repeat on per-sample basis
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
load([anstem '.globcov.' categdir '.N.mat'],'N','Ntot');
pat = load_struct([anstem '.patients.txt']);
P=[];
x = load_struct([anstem '.genomic.maf']);
for i=1:slength(pat), disp(i)
  xi = reorder_struct(x,strcmp(x.patient,pat.name{i}));
  P.title = pat.name{i};
  n = get_context_from_mafstruct(xi,[],categdir);
  Nn = [N{i}.callablebp(1:end-1) n];
  outstem = [anstem '.globcov.' categdir '.analysis.' pat.name{i}];
  [R H] = analyze_KaKs_data2(Nn,categdir,outstem,P);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  compare mutation rates in CpG islands, shores vs. sea
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
pat = load_struct([anstem '.patients.txt']);
load([anstem '.globcov.' categdir '.N.mat'],'N','Ntot');
x = load_struct([anstem '.genomic.maf']);
ntot = get_context_from_mafstruct(x,[],categdir);
Nn = [Ntot ntot];

% original (specific)  method
outstem = [anstem '.globcov.' categdir '.CpG_island.method1'];
R = analyze_CpG_island_data(Nn,categdir,outstem);

% new general method: allows user to arbitrarily configure:
%    
%     x  categories to exclude up-front (e.g. "bad", "any N")
%
%     q  (major division, shown as columns)
%     z  (minor division, shown as rows within each base-context category)
%     c  (base-context categories)
%         currently:  only two options
%                     1  (CpG/CG/AT)x(transition/transversion)
%                     2  melanoma categories
%         eventually: will accept a K struct from automatic mutation discovery

P=[];
P.exclude_grep = 'any N|bad';
P.major_title = 'CpG_zone';
P.major_grep = {'island','shore','sea','island|shore|sea'};
P.major_label = {'island','shore','sea','total'};
P.minor_title = 'transcript_zone';
P.minor_grep = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
P.minor_label = {'IGR','intron','UTR/pro','exon','total'};
outstem = [anstem '.globcov.' categdir '.CpG_island.method2'];
R = analyze_globcov_data(Nn,categdir,outstem,P);

% analyze cons vs. noncons,  transcript_zone
P=[];
P.exclude_grep = 'any N|bad';
P.major_title = 'regulatory_potential';
P.major_grep = {':cons','noncons','cons'};
P.major_label = {'RP','non-RP','total'};
P.minor_title = 'transcript_zone';
P.minor_grep = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
P.minor_label = {'IGR','intron','UTR/pro','exon','total'};
outstem = [anstem '.globcov.' categdir '.RP_vs_nonRP'];
R = analyze_globcov_data(Nn,categdir,outstem,P);
% per-sample
R = cell(slength(pat),1);
for i=1:slength(pat), disp(i)
  xi = reorder_struct(x,strcmp(x.patient,pat.name{i}));
  n = get_context_from_mafstruct(xi,[],categdir);
  Nn = [N{i}.callablebp(1:end-1) n];
  outstem = [anstem '.globcov.' categdir '.RP_vs_nonRP.' pat.name{i}];
  R{i} = analyze_globcov_data(Nn,categdir,outstem,P);
end
save([anstem '.globcov.' categdir '.RP_vs_nonRP.per_patient.analysis.mat'],'R');


% analyze good vs. bad,  transcript_zone
P=[];
P.exclude_grep = 'any N';
P.major_title = 'neighborhood';
P.major_grep = {'good','bad','good|bad'};
P.major_label = {'good','bad','total'};
P.minor_title = 'transcript_zone';
P.minor_grep = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
P.minor_label = {'IGR','intron','UTR/pro','exon','total'};
outstem = [anstem '.globcov.' categdir '.good_vs_bad'];
R = analyze_globcov_data(Nn,categdir,outstem,P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  also get stats using conservation29 (instead of regulatory3)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract global coverge from cbb
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
pat = load_struct([anstem '.patients.txt']);
dr = '/xchip/cga1/lawrence/crc/analysis/20100610/freeze';
pat.tumcbb = regexprep(pat.name,'(.*)',[dr '/$1/wgs/zonecov/$1-Tumor.covbb']);
pat.normcbb = regexprep(pat.tumcbb,'Tumor','Normal');
demand_file([pat.tumcbb;pat.normcbb]);
save_struct(pat,[anstem '.patients.txt']);
covdir = [anstem '.globcov.' categdir]; ensure_dir_exists(covdir);
outfile = regexprep(pat.name,'(.*)',[covdir '/$1.global_coverage.' categdir '.txt']);
P=[]; P.covbb_file_extension = 'txt';
get_global_coverage_stats([pat.tumcbb pat.normcbb],outfile,categdir,P);

% annotate genomic mutations by category
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
x = load_struct([anstem '.genomic.maf']);
x = setfield(x,categdir,get_context(x.chr,x.start,['/xchip/cga1/lawrence/db/' categdir]));
save_struct(x,[anstem '.genomic.maf']);

% process coverage stats
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
pat = load_struct([anstem '.patients.txt']);
covdir = [anstem '.globcov.' categdir];
N = cell(slength(pat),1);
for i=1:slength(pat), fprintf('%d/%d ',i,slength(pat));
  N{i} = load_struct_specify_string_cols([covdir '/' pat.name{i} '.global_coverage.' categdir '.txt'],2);
end,fprintf('\n');
Ntot = 0; for i=1:length(N), Ntot=Ntot+N{i}.callablebp(1:end-1); end
save([anstem '.globcov.' categdir '.N.mat'],'N','Ntot');

% analyze cons vs. noncons,  transcript_zone
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
load([anstem '.globcov.' categdir '.N.mat'],'N','Ntot');
x = load_struct([anstem '.genomic.maf']);
ntot = get_context_from_mafstruct(x,[],categdir);
Nn = [Ntot ntot];
P=[];
P.exclude_grep = 'any N|bad';
P.major_title = 'regulatory_potential';
P.major_grep = {':cons','noncons','cons'};
P.major_label = {'RP','non-RP','total'};
P.minor_title = 'transcript_zone';
P.minor_grep = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
P.minor_label = {'IGR','intron','UTR/pro','exon','total'};
outstem = [anstem '.globcov.' categdir '.RP_vs_nonRP'];
R = analyze_globcov_data(Nn,categdir,outstem,P);

% check that we get same result using all-in-memory method

% process coverage
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
W = loadCoverageFromWiggle([anstem '.samples_covered.wig']);
[Z C] = get_categs(categdir);
Ntot = calc_Ntot_from_wig_and_categ(W,Z,C);
save([anstem '.globcov.' categdir '.Ntot.from_all-in-memory_method.mat'],'Ntot');

% analyze cons vs. noncons,  transcript_zone
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
load([anstem '.globcov.' categdir '.Ntot.from_all-in-memory_method.mat'],'Ntot');
x = load_struct([anstem '.genomic.maf']);
ntot = get_context_from_mafstruct(x,[],categdir);
Nn = [Ntot ntot];
P=[];
P.exclude_grep = 'any N|bad';
P.major_title = 'regulatory_potential';
P.major_grep = {':cons','noncons','cons'};
P.major_label = {'RP','non-RP','total'};
P.minor_title = 'transcript_zone';
P.minor_grep = {'IGR','intron','UTR','exon','IGR|intron|UTR|exon'};
P.minor_label = {'IGR','intron','UTR/pro','exon','total'};
outstem = [anstem '.globcov.' categdir '.RP_vs_nonRP.using_all-in-memory_method'];
R = analyze_globcov_data(Nn,categdir,outstem,P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  coverage in CpG islands
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = load_struct('/xchip/cga1/lawrence/db/islands/CpG_islands.txt');
I.chr = convert_chr(I.chrom); I = reorder_struct(I,~isnan(I.chr));
I.start = str2double(I.chromStart); I.end = str2double(I.chromEnd);
I.len = I.end-I.start+1;
I.center = round((I.start+I.end)/2);
save('/xchip/cga1/lawrence/crc/analysis/20100610/crc9.CpG_islands.mat','I');

g = Genome(true); g.loadWiggle('/xchip/cga1/lawrence/crc/analysis/20100610/crc9.samples_covered.wig');
load('/xchip/cga1/lawrence/crc/analysis/20100610/crc9.CpG_islands.mat','I');

% CpG island vicinity coverage plot

mn = [200 275 350 450 650 900];
mx = [230 325 425 550 775 1100];
margin = 100000;
left = -margin; right = +margin;
C = zeros(length(mn),right-left+1);
for a=1:length(mn), disp(a), n=0; 
    Ia = reorder_struct(I,I.len>=mn(a) & I.len<=mx(a));
    for i=1:slength(Ia)
    chr = Ia.chr(i); st = Ia.center(i)+left; en = Ia.center(i)+right;
    try
      c = double(g.getContents(chr,st,en)); c(c==-128) = 0;
      if mean(c)<8 continue; end
      C(a,:) = C(a,:) + c'; n=n+1;
    catch me, end
  end, C(a,:) = C(a,:) / n; nn(a)=n;
end
figure(1);clf;
nsamps = 9; plot(C'/nsamps);
legend(str2cell(sprintf('%d-%d bp (n=%d)\n',[mn' mx' nn']')),'location','southeast');
dispmargin = 3000;
xlim([-dispmargin dispmargin]-left);x = get(gca,'xtick'); set(gca,'xticklabel',x+left);
xlabel('distance (bp) from center of CpG island','fontsize',20);
ylabel('fraction of samples callable','fontsize',20);
set(gca,'fontsize',20); set(gcf,'color',[1 1 1]);
text(1e5+2200,0.8,'size of CpG island','fontsize',20);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  mutation rates in introns of very long genes
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load mutations, coverage, genes, and categories (takes ~10 min)
tic
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
N = loadCoverageFromWiggle([anstem '.samples_covered.wig']);
G = get_refseq_introns_and_exons;
[Z C] = get_categs(categdir);
X = load_struct([anstem '.genomic.maf']);
X = make_numeric(X,{'chr','start','end'});
X = reorder_struct(X,strcmp('SNP',X.classification));
X.categnum = str2double(getfield(X,categdir));
toc

tic  % (takes ~1 hr)
[G IGR] = analyze_intron_rates(X,N,G,Z,C);
save([anstem '.intron_rates.with_versions.mat'],'G','IGR');
toc

anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'gc29izc65e';
load([anstem '.intron_rates.with_versions.mat'],'G','IGR');
P = []; P.legend = 'off'; P.group_xlabels = {'cons','noncons','total'};
B = plot_intron_rates(G,IGR,[5 8 2],P); title('all mutations');
plot_intron_rates(G,IGR,[6 9 3],P); title('CpG transitions only');
plot_intron_rates(G,IGR,[7 10 4],P); title('mutations other than CpG transitions');
plot_intron_rates(G,IGR,1);  % (for legend)







% what causes the effect:  long introns?  or long genomic footprints? (ans: both)

%load([anstem '.intron_rates.nothing_excluded.mat'],'G','B');
load([anstem '.intron_rates.exclude_bad.mat'],'G','B');

%mn = [0   2e3 4e3 1e4 2e4 4e4 1e5 2e5 4e5 6e5 8e5 1e6]';
%mx = [1e3 4e3 1e4 2e4 4e4 1e5 2e5 4e5 6e5 8e5 1e6 inf]';

mn = [0   1e4 2e4 4e4 1e5 2e5 5e5]';
mx = [1e4 2e4 4e4 1e5 2e5 5e5 inf]';


nbins = length(mn); Q = zeros(nbins); qn = zeros(nbins); qN = zeros(nbins);
for i=1:slength(G), if ~mod(i,1000), fprintf('%d/%d ', i, slength(G)); end
  fbn = find(G.footprint(i)<mx,1); if isempty(fbn), continue; end % __IGR__
  for j=1:G.n_introns(i)
    ibn = find(G.intron_ends{i}(j)-G.intron_starts{i}(j)+1<mx,1);
    Q(fbn,ibn) = Q(fbn,ibn) + 1;
    qN(fbn,ibn) = qN(fbn,ibn) + G.intron_N{i}(j);
    qn(fbn,ibn) = qn(fbn,ibn) + G.intron_n{i}(j);
end,end, fprintf('\n');
binname = str2cell(sprintf('%d-%d bp\n',[mn mx]'));

figure(1);clf;toshow = [4 5 6 7];
for r=1:length(toshow)
  subplot(1,length(toshow)+1,r);
  row = toshow(r); xn = qn(row,:); xN = qN(row,:); [rate ci] = binofit(xn,xN);
  barweb(1e6*rate',1e6*ci(:,1),1e6*ci(:,2),0.8,'intron size');
  ylabel('intron mutation rate / Mb','fontsize',20);
  title(['genes ' binname{row}],'fontsize',20);ylim([0 10]); set(gca,'fontsize',20);
end
legend(binname,'location','eastoutside');











P=[];
P.exclude_bad = true;
P.exclude_N = true;
[G B] = analyze_intron_rates(X,N,G,Z,C,P);
save([anstem '.intron_rates.exclude_bad.mat'],'G','B');

load([anstem '.intron_rates.exclude_bad.mat'],'G','B');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  look for local sequence determinants of mutations
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load mutations
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
x = load_struct([anstem '.genomic.maf']);
outname = [anstem '.mutation_vicinity_analysis.txt'];
analyze_mutation_vicinity(x,outname);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  look for interesting rearrangement that have breakpoints in exons
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
pat = load_struct([anstem '.patients.txt']);
x = cell(slength(pat),1);
for i=1:slength(pat), disp(i)
  tmp = load(['/xchip/cga1/lawrence/crc/analysis/20100610/freeze/' pat.name{i} ...
     '/wgs/ra/' pat.name{i} '.dRanger_results.detail.all.mat']);
  x{i} = tmp.X;
end
x = concat_structs(x);
xall=x;
x = reorder_struct(xall,xall.somatic_score>0);
x = dRanger_annotate_sites(x);  % fix IGR annotations
save_struct(x,'/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic_rearrs_with_fixed_annotations.txt');
save('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic_rearrs_with_fixed_annotations.mat','x');

x = reorder_struct(xall,xall.somatic_score>=3);
x = dRanger_annotate_sites(x);  % fix IGR annotations
save_struct(x,'/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic3_rearrs_with_fixed_annotations.txt');
save('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic3_rearrs_with_fixed_annotations.mat','x');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  look for patterns of microhomology vs. foreign sequence at RA bkpts
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic_rearrs_with_fixed_annotations.mat','x');
x = reorder_struct(x,x.BPresult==1);
x = reorder_struct(x,x.somatic_score>=3);
xcount(x.individual,x.T_lenhomology)
xcount(x.individual,x.T_lenforeign)

x = reorder_struct(x,x.somatic_score>=8);
xcount(x.individual,x.T_lenhomology)


load('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic_rearrs_with_fixed_annotations.mat','x');
x = reorder_struct(x,x.somatic_score>=3);
x.type1=regexprep(x.site1,'^(Exon|Intron|Promoter|5''-UTR|3''-UTR|IGR).*','$1')
x.type2=regexprep(x.site2,'^(Exon|Intron|Promoter|5''-UTR|3''-UTR|IGR).*','$1')
xcount(x.type1,x.type2)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  look for recurrent breakpoints
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic_rearrs_with_fixed_annotations.mat','x');
x = reorder_struct(x,x.somatic_score>=3); nx = slength(x);
x.ino = str2double(regexprep(x.individual,'CRC-(\d*)','$1'));
b = [x.chr1 x.pos1 x.str1 x.ino (1:nx)' ones(nx,1);x.chr2 x.pos2 x.str2 x.ino (1:nx)' 2*ones(nx,1)];
b = sortrows(b);
d = diff(b(:,2)); d(b(1:end-1,1)~=b(2:end,1)) = inf;
d(b(1:end-1,4)==b(2:end,4)) = nan;
b = [[d;nan] b];
idx = find(b(:,1)>=1e4 & b(:,1)<=3e4)
idx2 = sort([idx;idx+1]);
b(idx2,:)

% permutations
for i=1:10000
  q = 3e9*rand(1368,1);
  dq = diff(sort(q));
  s(i) = sum(dq<1e4);
end

% genes hit in multiple samples

y = x;
exclude_igr = true;
if exclude_igr
  idx=grep('^IGR',y.site1,1); y.gene1(idx) = repmat({'---'},length(idx),1);
  idx=grep('^IGR',y.site2,1); y.gene2(idx) = repmat({'---'},length(idx),1);
end
[u ui uj] = unique([y.gene1;y.gene2]);
q = false(length(u),9);
for i=1:9, disp(i)
  yi = reorder_struct(y,y.ino==i+1);
  for j=1:length(u)
    q(j,i) = any(strcmp(yi.gene1,u{j})) | any(strcmp(yi.gene2,u{j}));
end,end
sq = sum(q,2);
histc(sq,1:9)
%G = load_gene_lengths;idx = listmap(u,G.name);

v = false(slength(y));
threshold = 100000;
for i=1:slength(y), if ~mod(i,100), disp(i), end, for j=1:slength(y)
  v(i,j)=(y.ino(i)~=y.ino(j) & y.chr1(i)==y.chr1(j) & y.chr2(i)==y.chr2(j) &...
    abs(y.pos1(i)-y.pos1(j))<threshold & abs(y.pos2(i)-y.pos2(j))<threshold);
end,end
[a1,a2,a3]=find(v);
unique([y.gene1(a1);y.gene2(a1)])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  combine mutations, indels, and rearrangements
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% point mutations and indels
P = [];
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
P.mutfile = [anstem '.with_null_category.maf'];
P.covfile = [anstem '.coverage.C1.mat'];
P.catfile = [anstem '.mutcategs.with_null_category.txt'];
P.patlist = [anstem '.patients.txt'];
P.genelist = '/xchip/tcga_scratch/lawrence/capture/weRefseq_genelist.txt';
M = load_all_mutation_data2(P);
M = analyze_mutation_data(M,anstem,P);

% rearrangements
R = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic_rearrs_with_fixed_annotations.txt');
R.type1=regexprep(R.site1,'^(Exon|Intron|Promoter|5''-UTR|3''-UTR|IGR).*','$1')
R.type2=regexprep(R.site2,'^(Exon|Intron|Promoter|5''-UTR|3''-UTR|IGR).*','$1')
Rg = reorder_struct(R,union(grepv('IGR',R.type1,1),grepv('IGR',R.type2,1)));
map1 = listmap(Rg.gene1,M.gene.name); map2 = listmap(Rg.gene2,M.gene.name);
M.n_rearrs = zeros(M.ng,1); M.n_rearr_samps = zeros(M.ng,1);
for i=1:M.ng
  idx = find(map1==i | map2==i);
  M.n_rearrs(i) = length(idx);
  M.n_rearr_samps(i) = length(unique(R.individual(idx)));
end

% list most rearranged genes and their nums of mutations
[tmp ord] = sort(M.n_rearrs,'descend');
if 1, fprintf('%-15s\t%s\t%s\t%s\t%s\t%s\n','gene','nrearrs','nrsamps','nmuts','p','q');
for i=1:30, j=ord(i);
  fprintf('%-15s\t%d\t%d\t%d\t%0.0d\t%0.0d\n',M.gene.name{j},M.n_rearrs(j),M.n_rearr_samps(j),...
     sum(M.n_used(j,:)),M.Prob(j),M.Q(j));
end, end

% list most mutated genes and their num of rearrs
[tmp ord] = sort(M.Prob);
if 1, fprintf('%-15s\t%s\t%s\t%s\t%s\t%s\n','gene','nrearrs','nrsamps','nmuts','p','q');
for i=1:30, j=ord(i);
  fprintf('%-15s\t%d\t%d\t%d\t%0.0d\t%0.0d\n',M.gene.name{j},M.n_rearrs(j),M.n_rearr_samps(j),...
     sum(M.n_used(j,:)),M.Prob(j),M.Q(j));
end, end

% look at genes with both mutations and rearrangements
idx = find(M.n_rearrs>0 & sum(M.n_used,2)>0);
[tmp ord] = sort(M.Prob(idx));
if 1, fprintf('%-15s\t%s\t%s\t%s\t%s\t%s\n','gene','nrearrs','nrsamps','nmuts','p','q');
for i=1:30, j=idx(ord(i));
  fprintf('%-15s\t%d\t%d\t%d\t%0.0d\t%0.0d\n',M.gene.name{j},M.n_rearrs(j),M.n_rearr_samps(j),...
     sum(M.n_used(j,:)),M.Prob(j),M.Q(j));
end, end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 %%%%
%%%%  TABLE 1        %%%%
%%%%                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%   clinical info (diagnosis, age, gender, tumor%)
%   purity and ploidy estimates from Scott
T = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/table_from_adam.txt',[],[],true);
T = rename_field(T,{'firehoseid','tumor'},{'individual','pathologist_purity_estimate_pct'});
T = order_fields_first(T,{'individual'});
S = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/purity_ploidy_from_scott.txt',[],[],true);
S = make_numeric(S,{'purity','ploidy'});
T.computationally_estimated_purity_pct = 100*mapacross(T.collaboratortumorsample,S.sample_id,S.purity);
T.computationally_estimated_ploidy = mapacross(T.collaboratortumorsample,S.sample_id,S.ploidy);
T = reorder_struct(T,ismember(T.individual,str2cell(sprintf('CRC-%04d\n',2:10))));

% global genome coverage
for i=1:slength(T)
  z = load_struct(['/xchip/cga1/lawrence/crc/analysis/20100610/freeze/'...
    T.individual{i} '/wgs/zonecov/' T.individual{i} '_global_coverage_by_zone.txt']);    
  z = make_numeric(z,{'tdepth','ndepth','tseqbp','nseqbp','callablepct'});
  T.tumor_sequenced_bp(i,1) = z.tseqbp(end);
  T.normal_sequence_bp(i,1) = z.nseqbp(end);
  T.tumor_num_reads(i,1) = round(z.tseqbp(end)/101);
  T.normal_num_reads(i,1) = round(z.nseqbp(end)/101);
  T.tumor_coverage_depth(i,1) = z.tdepth(end);
  T.normal_coverage_depth(i,1) = z.ndepth(end);
  T.genome_callable_pct(i,1) = z.callablepct(end);
end

% number of genomic mutations (and rate)
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
x = load_struct([anstem '.genomic.maf']);
for i=1:slength(T), T.num_genomic_mutations(i,1) = sum(strcmp(x.patient,T.individual{i})); end
T.mutation_rate_whole_genome_per_mb = 1e6 * T.num_genomic_mutations ./ (3e9*T.genome_callable_pct/100);

% number of coding mutations (nonsilent and silent and indels)
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
x = load_struct([anstem '.maf']);
for i=1:slength(T)
  idx = find(strcmp(x.patient,T.individual{i}));
  T.num_coding_mutations(i,1) = length(idx);
  T.num_coding_mutations_point_nonsilent(i,1) = length(grep('Mis|Nons|Read-through|Splice',x.type(idx),1));
  T.num_coding_mutations_point_silent(i,1) = length(grep('Syn',x.type(idx),1));
  T.num_coding_mutations_indel_inframe(i,1) = length(grep('In_Frame',x.type(idx),1));
  T.num_coding_mutations_indel_frameshift(i,1) = length(grep('Frame_Shift',x.type(idx),1));
end

% number of rearrangements
R = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/freeze/ra/CRC_somatic_rearrs_with_fixed_annotations.txt');
R = make_numeric(R,{'chr1','chr2','somatic_score'});
R = reorder_struct(R,R.somatic_score>=3);
for i=1:slength(T)
  idx = find(strcmp(R.individual,T.individual{i}));
  T.num_rearrangements(i,1) = length(idx);
  T.num_rearrangements_intrachomosomal(i,1) = sum(R.chr1(idx)==R.chr2(idx));
end
T.num_rearrangements_interchomosomal = T.num_rearrangements - T.num_rearrangements_intrachomosomal;

save('/xchip/cga1/lawrence/crc/analysis/20100610/Table1.mat','T');
save_struct(T,'/xchip/cga1/lawrence/crc/analysis/20100610/Table1.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 %%%%
%%%%  FIGURE 3       %%%%
%%%%                 %%%%
%%%%  TCF7L2/VTI1A   %%%%
%%%%  rearrangement  %%%%
%%%%                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% (done mostly in powerpoint)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  %%%%
%%%%  FIGURE 1        %%%%
%%%%                  %%%%
%%%%     mutation     %%%%
%%%% rates and biases %%%%
%%%%                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIGURE 1a
% rates per sample, broken down by base category
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
pat = load_struct([anstem '.patients.txt']);
load([anstem '.globcov.' categdir '.RP_vs_nonRP.per_patient.analysis.mat'],'R');
[tmp ui] = unique(R{1}.base_context); bc = R{1}.base_context(sort(ui));
rate = nan(slength(pat),length(bc),3);
for i=1:slength(pat)
  r = reorder_struct(R{i},strcmp('total',R{i}.regulatory_potential) & strcmp('total',R{i}.transcript_zone));
  for j=1:length(bc)
    idx = find(strcmp(bc{j},r.base_context));
    rate(i,j,1) = r.mutrate(idx);
    rate(i,j,2) = r.ci_low(idx);
    rate(i,j,3) = r.ci_high(idx);
end,end
% sort by total rate
[tmp ord] = sort(rate(:,end,1));
rate = rate(ord,:,:);
patname = pat.name(ord);
bc = regexprep(bc,'_',' ');
% draw figure
figure(1),clf, fs = 22;
subplot('position',[0.1 0.1 0.12 0.8]);
barweb(rate(:,1,1)',rate(:,1,3)',rate(:,1,2)',0.8,bc(1)); ylabel('mutations per million sites','fontsize',fs);
set(gca,'fontsize',fs);
subplot('position',[0.26 0.1 0.67 0.8]);
barweb(rate(:,2:end,1)',rate(:,2:end,3)',rate(:,2:end,2)',0.8,bc(2:end)); ylabel('mutations per million sites','fontsize',fs);
set(gca,'fontsize',fs);
set(gcf,'position',[1000 20 2399 904],'color',[1 1 1]);
legend(patname,'location','northwest');

%%%% 2010-07-27
%%%% Alternate version, with categories stacked, and denominator = total bases sequenced
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
pat = load_struct([anstem '.patients.txt']);
load([anstem '.globcov.' categdir '.RP_vs_nonRP.per_patient.analysis.mat'],'R');
[tmp ui] = unique(R{1}.base_context); bc = R{1}.base_context(sort(ui));
N = nan(slength(pat),1);
n = nan(slength(pat),length(bc),1);
for i=1:slength(pat)
  r = reorder_struct(R{i},strcmp('total',R{i}.regulatory_potential) & strcmp('total',R{i}.transcript_zone));
  N(i) = r.Ncov(strcmp('total',r.base_context));
  for j=1:length(bc)
    idx = find(strcmp(bc{j},r.base_context));
    n(i,j) = r.n_muts(idx);
end,end
rate = bsxfun(@rdivide,n,N);
% sort by total rate
[tmp ord] = sort(rate(:,end));
rate = rate(ord,:);
patname = pat.name(ord);
% draw figure
figure(1),clf
bar(1e6*rate(:,1:end-1),'stacked')
ylabel('mutations per megabase','fontsize',20);
set(gca,'xticklabel',patname,'fontsize',20);
set(gcf,'position',[1000 20 1900 904],'color',[1 1 1]);
legend(bc(1:end-1),'location','northwest','interpreter','none');





% FIGURE 1b
% IGR vs. intron vs. UTR/promoter vs. exon
%     in total, RP, non-RP
%     and in all mutations, CpG transitions, and mutations other than CpG transitions

categdir = 'grizc65e'; ww = 'RP';      % 3-mammal alignment ("RP vs. non-RP")
categdir = 'gc29izc65e'; ww = 'cons';  % 29-mammal alignment ("cons vs. non-cons")
% load data
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
load([anstem '.globcov.' categdir '.RP_vs_nonRP.analysis.mat'],'R');
P=[]; P.cons_term = ww;
draw_cons_vs_noncons_plot(R,P);

% FIGURE 1c
% IGR vs. intron vs. UTR/promoter vs. exon
%     in CpG_island, CpG_shore, CpG_sea, total
%     and in all mutations, CpG transitions, and mutations other than CpG transitions

% load data
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9'; categdir = 'grizc65e';
load([anstem '.globcov.' categdir '.CpG_island.method2.analysis.mat'],'R');
P=[]; P.fontsize = 12;
draw_CpG_island_rate_plot(R,P);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  combined analysis with Hopkins discovery set (Wood et al.)
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load Wood et al. target list (has been LiftOver'ed from hg17 to hg18)
T = load_struct('/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS1_liftover.txt');
T = reorder_struct(T,grepv('unmapped',T.hg18,1));
tmp = parse(T.hg18,'(chr.*):(\d*)-(\d*)',{'chr','start','end'},2:3);
T = merge_structs({T,tmp});
T.chr = convert_chr(T.chr);

% make wiggle file of targeted bases
g = Genome;
flank = 2;
for i=1:slength(T), if ~mod(i,10000), fprintf('%d/%d ',i,slength(T)); end
  g.setContents(T.chr(i),T.start(i)-flank,T.end(i)+flank,1);
end, fprintf('\n');
fprintf('Total territory: %ld\n', g.getTotContents());
g.saveWiggle('/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS1_hg18.wig');

% convert to coverage file
targlist = ['/xchip/tcga_scratch/lawrence/capture/' ...
           'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP_noheader.txt'];
categfile = '/xchip/tcga_scratch/lawrence/db/context65/weRefseq_terr_only.wig';
mincateg = 1; maxcateg = 65;
names = {'Vog'};
wigfiles = {'/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS1_hg18.wig'};
outfiles = {'/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS1_hg18.weRefseq.somatic_coverage.txt'};
extract_from_wig(names,targlist,outfiles,categfile,mincateg,maxcateg,wigfiles);

% check on coverage plot (to make sure LiftOver worked)
C = []; C.sample = []; C.sample.name = {'VOG'};
C.sample.covfile = {'/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS1_hg18.weRefseq.somatic_coverage.txt'};
C.ns = slength(C.sample);
C.file.targ = ['/xchip/tcga_scratch/lawrence/capture/' ...
           'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP_noheader.txt'];
C.file.categdir ='/xchip/tcga_scratch/lawrence/db/context65'; C.ncat = 65;
C = load_coverage(C);
K = load_struct('/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/categs_CpG_CG_AT_Indel.txt');
C1 = collapse_coverage_categories(C,K);
C1 = compute_fraction_coverage(C1);
save('/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS1_hg18.C1.mat','C1','-v7.3');
capture_covplot(C1);   % ---> looks good!

% combined patient list
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9';
pat = load_struct([anstem '.patients.txt']);
vogpat = load_struct('/xchip/tcga/gbm/analysis/lawrence/vog/brco/11colon_discovery_patients.txt')
vogpat.covfile = repmat({'/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS1_hg18.weRefseq.somatic_coverage.txt'}, ...
                        slength(vogpat),1);
pat = concat_structs_keep_all_fields({pat,vogpat});
save_struct(pat,'/xchip/cga1/lawrence/crc/analysis/20100610/crc9_plus_vog11.patients.txt');

% make combined coverage struct and collapse categories
anstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9_plus_vog11';
pat = load_struct([anstem '.patients.txt']);
C = []; C.sample = pat; C.ns = slength(C.sample);
C.file.targ = ['/xchip/tcga_scratch/lawrence/capture/' ...
           'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP_noheader.txt'];
C.file.categdir ='/xchip/tcga_scratch/lawrence/db/context65'; C.ncat = 65;
C = load_coverage(C);
K = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/crc9.mutcategs.txt');
C1 = collapse_coverage_categories(C,K);
C1 = compute_fraction_coverage(C1);
save([anstem '.coverage.C1.mat'],'C1','-v7.3');
save([anstem '.coverage.C.mat'],'C','-v7.3');

% mutations: make into a MAF
M = load_struct('/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS3_liftover.txt');
M = reorder_struct(M,strcmp(M.screen,'Discovery') & strcmp(M.tumor_type,'Colorectal'));
length(unique(M.tumor)) % 11
M = make_numeric(M,{'start','end'});
tmp = parse(M.hg18,'\d*(([ACGT]>[ACGT])|((del|ins|dup)[ACGT]*))(( \(homozygous\))?)',{'change','hom'});
M = merge_structs({M,tmp});
M.ishom = ~cellfun('isempty',M.hom);
Mi = reorder_struct(M,grepi('indel',M.type,1));
tmp = parse(Mi.change,'(dup|del|ins)([ACGT]*)',{'kind','seq'});
Mi = merge_structs({Mi,tmp});
Mdel = reorder_struct(Mi,strcmp('del',Mi.kind));
% fix strandedness for deletions
for i=1:slength(Mdel),Mdel.r{i,1}=genome_region(Mdel.chr{i},Mdel.start(i),Mdel.end(i));end  % get actual reference
idx = find(~strcmpi(Mdel.r,Mdel.seq));
Mdel.seq(idx) = rc(Mdel.seq(idx));
Mdel.ref = Mdel.seq;
Mdel.tum1 = repmat({'-'},slength(Mdel),1);
Mdup = reorder_struct(Mi,strcmp('dup',Mi.kind));
pos = Mdup.end; Mdup.start = pos; Mdup.end = pos+1;
Mdup.ref = repmat({'-'},slength(Mdup),1);
Mdup.tum1 = Mdup.seq;
% hard to fix strandedness for duplications / insertions
Mins = reorder_struct(Mi,strcmp('ins',Mi.kind));
Mins.ref = repmat({'-'},slength(Mins),1);
Mins.tum1 = Mins.seq;
Mp = reorder_struct(M,grepvi('indel',M.type,1));
tmp = parse(Mp.change,'([ACGT])>([ACGT])',{'ref','alt'});
Mp = merge_structs({Mp,tmp});
% fix strandedness for point mutations
for i=1:slength(Mp),Mp.r{i,1}=genome_region(Mp.chr{i},Mp.start(i));end  % get actual reference
idx = find(~strcmpi(Mp.r,Mp.ref));
Mp.ref(idx) = rc(Mp.ref(idx)); Mp.alt(idx) = rc(Mp.alt(idx));
Mp.tum1 = Mp.alt;
M2 = concat_structs_keep_all_fields({Mdel,Mins,Mdup,Mp});
M2.tum2 = M2.ref; M2.tum2(M2.ishom) = M2.alt(M2.ishom);
x = []; x.col1 = repmat({'36'},slength(M2),1);
x.col2 = M2.chr; x.col3 = M2.start; x.col4 = M2.end;
x.col5 = M2.ref; x.col6 = M2.tum1; x.col7 = M2.tum2;
x.col8 = M2.tumor; x.col9 = regexprep(M2.tumor,'(.*)','$1-normal');
save_struct_noheader(x,'/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS3_liftover.maflite');
annotate_maflite('/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS3_liftover.maflite',...
   '/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS3_liftover.maf.annotated','hg18');

% combine Vog mutations with our data
x = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/crc9.maf');
y = load_struct('/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS3_liftover.maf.annotated');
y.Hugo_Symbol = y.gene; % (to avoid empty first column: causes problems for load_struct)
y.patient = y.tumor_barcode;
y.dataset = repmat({'VOG'},slength(y),1);
z = concat_structs_keep_all_fields({x,y});
z = reorder_struct(z,grepv('UTR|Intron|Promoter|miRNA|IGR|Non-coding|Non-mutation',z.type,1));
z.context = get_context(z.chr,z.start,'/xchip/tcga_scratch/lawrence/db/context65');
z = collapse_adjacent_mutations(z);  % (did nothing)
save_struct(z,'/xchip/cga1/lawrence/crc/analysis/20100610/crc9_plus_vog11.maf');



%%%%%%%%%%%%%%%%%%%
% standard MutSig %
%%%%%%%%%%%%%%%%%%%

P=[];
P.covfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.coverage.C1.mat';
P.mutfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.with_null_category.maf';
P.patlist = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.patients.txt';
P.catfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.mutcategs.with_null_category.txt';
P.genelist = '/xchip/tcga_scratch/lawrence/capture/weRefseq_genelist.txt';
outstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.run2';
isetname = 'crc9.run2';
perform_mutsig_analysis(isetname,outstem,P)


%%%%%%%%%%%%%%%%%%%
% combined MutSig %
%%%%%%%%%%%%%%%%%%%

P=[];
P.covfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9_plus_vog11.coverage.C1.mat';
P.mutfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9_plus_vog11.maf';
P.patlist = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9_plus_vog11.patients.txt';
P.catfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.mutcategs.with_null_category.txt';
P.genelist = '/xchip/tcga_scratch/lawrence/capture/weRefseq_genelist.txt';
outstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9_plus_vog11';
isetname = 'crc9_plus_vog11';
perform_mutsig_analysis(isetname,outstem,P)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard MutSig on genes that made it to Hopkins validation phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get list of all genes with at least one non-synonymous mutation in the Hopkins discovery set

y = load_struct('/xchip/tcga/gbm/analysis/lawrence/vog/brco/tableS3_liftover.maf.annotated');
y = reorder_struct(y,grep('Del|Ins|Missense|Nonsense|Read-through',y.type,1));
g =[]; g.name = unique(y.gene);  % 729 genes
save_struct(g,'/xchip/tcga/gbm/analysis/lawrence/vog/brco/colon_validation_genelist.txt');

P=[];
P.covfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.coverage.C1.mat';
P.mutfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.with_null_category.maf';
P.patlist = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.patients.txt';
P.catfile = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.mutcategs.with_null_category.txt';
P.genelist = '/xchip/tcga/gbm/analysis/lawrence/vog/brco/colon_validation_genelist.txt';
outstem = '/xchip/cga1/lawrence/crc/analysis/20100610/crc9.hopvalgenes';
isetname = 'crc9.hopvalgenes';
perform_mutsig_analysis(isetname,outstem,P)




%%%%% recurrent mutations for manual review

m = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/crc9.maf');
x = load_struct('/xchip/cga1/lawrence/crc/analysis/20100610/RECURRENT MUTATIONS for man review.txt');

g = unique(x.Hugo_Symbol);
m2 = reorder_struct(m,ismember(m.Hugo_Symbol,g));
m2 = reorder_struct(m2,grepv('Silent|Syn',m2.Variant_Classification,1));

mi = reorder_struct(m2,grepi('In|Del',m2.Variant_Classification,1));
save_struct(mi,'/xchip/cga1/lawrence/crc/analysis/20100610/indels_for_review.txt');



