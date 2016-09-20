function [XY,X]=plotMutationSpectrumCategLego1(MAF,COV,PNGFILE,P)
%function [XY]=plotMutationSpectrumCategLego1(MAF,COV,PNGFILE,P)
%
% plot and save mutation spectrum lego plot for one sample based on maf file
% and coverasge files in standard firehose location
%
% inputs:   MAF: filename or MAF structure with fields:
%               .ref_context
%               .Tumor_Sample_Barcode
%               .Matched_Norm_Sample_Barcode
%               .Variant_Type (='SNP' only, DNPs or Indels are not included in the plot)
%               .Reference_Allele
%               .Tumor_Seq_Allele1
%           COV: Coverage file from mutsig, 'exome' for precomputed, 'unit' for no coverage
%           PNGFILE: output lego plot png file
%           P: optional parameters
%              .zscale = true (print zscale on main rate plot)
%              .label = label for lego to override the first Tumor_Sample_Barcode
%              
%
% outputs:  XY: structure for lego plot
%               XY.n:  8x12 matrix of mutation counts per lego bin
%               XY.ncov: 8x12 matrix of base coverage per lego bin
%               XY.cat:  8x12 cell matrix of context+mutation categories per lego bin
%					     eg. G-C>A-T is a C>A mutation after a G before a T
%               XY.col:  8x12x3 matrix of RGB colors per lego bin'
%
% CS:  19Sep2012, based on plotMutationSpectrumCategLegos

fprintf('plotMutationSpectrumCategLegos\n');

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'zscale',true);
P = impose_default_value(P,'subplot',false);

MAFFILE='';
if ~isstruct(MAF) && exist(MAF,'file')
  fprintf('MAF file:\t %s\n',MAF);
  MAFFILE=MAF;
  MAF=load_table(MAF);
else
  fprintf('MAF structure length:\t %d\n',length(MAF.Variant_Type));
  k=find(ismember(MAF.Variant_Type,'SNP'));
  if (length(k)<length(MAF.Variant_Type))
      fprintf('MAF structure SNVs:\t %d\n',length(k));
      MAF=trimStruct(MAF,k);
  end
end
P = impose_default_value(P,'label',MAF.Tumor_Sample_Barcode{1});

% categ struct
base='ACGT';
n=0; C=[];
categs.name=repmat({''},4*4*4*13+13,1);
for i=1:4; for j=1:4;  for k=1:4
   n=n+1;
   C.name(n,1)={[base(i) ' in ' base(j) '_' base(k)]};
   C.id(n,1)=n;
end;    end ;end
n=n+1;
C.name(n,1)={'any N'};
C.id(n,1)=n;
C.N=length(C.id);

[XY,X,C1]=MutationSpectrumCategs(MAF,COV);
if exist(COV,'file')
  fprintf('COV file:\t %s\n',COV);
else
  fprintf('use precomputed coverage:\t %s\n',COV);
end      
fprintf('write png image to:\t %s\n',PNGFILE);


%%
% artifact context category
ART='C-C>A-G';
kART=find(ismember(XY.cat,{ART}));
nART=XY.n(kART);
fART=XY.n(kART)/sum(XY.n(:));
ARTX='C-C>A-';
ARTX='-C>A-';  % exclude all C>A for non-artifact selection
kARTX=strfindk(XY.cat,ARTX,'v');
CpGT='C>T-G';
kCpGT=strfindk(XY.cat,CpGT);

% rate per Mb per sample 
Np=1
z=1e6;

if (round(Np*prod(size(XY.ncov)))==round(sum(XY.ncov(:))))
    z=1;
    zlab=' SNV count  ';
    fprintf('mutation count lego \n');
else    
    zlab=' SNV rate/Mb ';
    fprintf('mutation rate lego\n');

end
%if (~P.visible), figure('Visible','Off'); end
if (~P.subplot)
    clf
    fp=get(gcf,'position');
    set(gcf,'position',[fp(1:2) 924 604])
end

bar3_with_colors(z*XY.n./XY.ncov,XY.col)
if (~P.zscale)
	set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);
else
	set(gca,'xticklabel',[],'yticklabel',[],'ylim',[0.5 8.5]);
	zlabel(zlab);	
end
set(gcf,'color',[1 1 1])
h0=gca;
xp=0.05;
if (P.zscale), xp=0.125; end
text(xp,0.95,[P.label ' '],'units','normalized','fontsize',11,'interpreter','none');
text(xp,0.9,sprintf('n=%d',X.N),'units','normalized','fontsize',11)    

q=regexp(X.CT,'-','split');
q=cellfun(@(x) x(2),q);
SNVT={'C>T';'C>A';'C>G';'A>G';'A>C';'A>T'};
nc=zeros(6,1);
for i=1:length(SNVT)
  k=strfindk(q(:),SNVT{i});
  nc(i)=length(k)+1e-5;
end
yp=.6; 
if (P.zscale), yp=0.7; end
axes('Position', [0.05,yp, .2, .2]);
h=pie(nc,repmat({''},length(nc),1));
colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7;];

for i=1:6
  set(h(2*i-1),'facecolor',colors(i,:))
end;
h1=legend(regexprep(SNVT,'>','>'));
set(h1,'Position',[.25 .6 .1 .2])
q=regexprep(XY.cat(1:4,1:4),'-','_'); 
box1=cellfun(@(x) [x([1 2 7]) ' '],q,'UniformOutput',false);
if (~P.zscale)
	axes('Position', [0.53,0.1, .5, .5]);
else
	axes('Position', [0.57,0.05, .5, .5]);
end
h=bar3(ones(4,4),'w');
set(h,'EdgeAlpha',0.25);
set(gca,'zlim',[0 50]);
%bar3([0 1],[0 1],[0 1],'color','w')
for i=1:4,for j=1:4
  hh(i,j)=text((j-0.35),(i)+0.15,4,box1(i,j),'horizontalAlign','left','EraseMode','none','interpreter','none'); %,'units','normalized')
end,end
set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);
set(h1,'OuterPosition',[0.7871 0.5546 0.1616 0.2670])
set(h1,'Position', [0.7942 0.5594 0.1509 0.2622])

%%
if ~exist('PNGFILE','var'), 
    PNGFILE=[ '~/mutation_context_lego' P.label]; 
end
saveas(gcf,[PNGFILE '.png'],'png')


%%
function test

MAF='/local/cga-fh/cga/Sigma_Pediatric/Individual_Set/PR_SIGMA_Rhabdoid_WGS/Individual/Rhabdoid-07-221/jobs/wgs/mut/annotated/Rhabdoid-07-221-Tumor.maf.annotated'
COV='unit'
PNGFILE='lego.count'
[XY]=plotMutationSpectrumCategLego1(MAF,COV,PNGFILE)

MAF='/local/cga-fh/cga/Sigma_Pediatric/Individual_Set/PR_SIGMA_Rhabdoid_WGS/Individual/Rhabdoid-09-046A/jobs/capture/mut/annotated/Rhabdoid-09-046A-Tumor.maf.annotated'
COV='exome'
PNGFILE='lego.rate'
[XY]=plotMutationSpectrumCategLego1(MAF,COV,PNGFILE)

MAF='/local/cga-fh/cga/Sigma_Pediatric/Individual_Set/PR_SIGMA_Rhabdoid_Capture/jobs/mutsig1.5_oxog/PR_SIGMA_Rhabdoid_Capture.final_analysis_set.maf'
COV='/local/cga-fh/cga/Sigma_Pediatric/Individual_Set/PR_SIGMA_Rhabdoid_Capture/jobs/mutsig1.5_oxog/PR_SIGMA_Rhabdoid_Capture.coverage.mat'
PNGFILE='~/Desktop/lego'
[XY]=plotMutationSpectrumCategLego1(MAF,COV,PNGFILE)
