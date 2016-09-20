function [XY,X,C1]=plotMutationSpectrumCategLegos(WORKSPACE,SET,MUTSIG,OUT,P)
%function [XY,X,C1]=plotMutationSpectrumCategLegos(WORKSPACE,SET,MUTSIG,OUT,P)
%function [XY,X,C1]=plotMutationSpectrumCategLegos(MAFFILE,OUT,P)
%
% plot and save mutation spectrum lego figures based on MutSig maf file
% and coverasge files in standard firehose location
%
% inputs:   WORKSPACE: workspace name
%					   (under xchip/cga1/firehose_output)
%           SET: workspace name
%					   (under /xchip/cga1/firehose_output/<WORKSPACE>/Individual_Set)
%           MUTSIG: Mutsig version (mutsig1.5 default)
%					   (under /xchip/cga1/firehose_output/<WORKSPACE>/Individual_Set/<SET>)
%           OUT: output area  (~/<WORKSPACE>/plots/<datestr(now, 'ddmmmyyyy')> default
%           P: parameter control structure
%              X = maf structure (the whole structure)
%              C1 = mutsig coverage  structure (needs only C1.orig_cov)
%                   can also be the string "exome" or "genome", in which case a default coverage distribution will be used,
%                   based on a precomputed average from many previous well-covered samples.
%                   The number of patients will be estimated from length(unique(maf.patient))--
%                     this only matters for the absolute rate estimate, displayed if zscale==true.
%                   A custom coverage file with a vector of 65 values can override 'genome' or 'exome'
%                   
%
%              zscale = true (print zscale on main rate plot)
%              printrates = true (print rate rows in text output)
%              FH_AREA = optional input area for maf files and coverage mat files
%              unix = skip any unix commands (hangs on some enviroments)
%
% alternative input structure:   
%
%           MAFFILE: direct path to MAF file
%           OUT: as above
%           P: as above
%
% outputs:  XY: structure for lego plot
%               XY.n:  8x12 matrix of mutation counts per lego bin
%               XY.ncov: 8x12 matrix of base coverage per lego bin
%               XY.cat:  8x12 cell matrix of context+mutation categories per lego bin
%					     eg. G-C>A-T is a C>A mutation after a G before a T
%               XY.col:  8x12x3 matrix of RGB colors per lego bin
%           X: MAF struction with Mutsig fields & added CT field with mutation category cellstr
%           C1:coverge structure stripped down to a few fields including orig_cov
%
% CS:  04Mar2012, based on Mike Lawrence's draw_mutation_spectrum_3d_barplot
% ML:  15Mar2012, added option to use precomputed average coverage ("exome" or "genome")
% ML:  08May2012, added alternative parameter structure, to point directly to MAF

fprintf('plotMutationSpectrumCategLegos\n');

if exist(WORKSPACE,'file')
    if (nargin>2)
        P=MUTSIG;
    end
    OUT=SET;
end
if isstruct(MUTSIG)
    P=MUTSIG;
    MUTSIG='';
end
    

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'zscale',false);
P = impose_default_value(P,'printrates',false);
P = impose_default_value(P,'unix',true);
P = impose_default_value(P,'FH_AREA','/local/cga-fh/cga');
% P = impose_default_value(P,'visible',false);


if exist(WORKSPACE,'file')
  % alternative parameter structure: plotMutationSpectrumCategLegos(MAFFILE,OUT,P)
  maffile = WORKSPACE;
  covfile = 'exome';
  clear WORKSPACE;
  if exist('SET','var')
    OUT = SET;
    clear SET;
  end
  if exist('MUTSIG','var')
    %P = MUTSIG;
    clear MUTSIG;
    if isfield(P,'C1'), covfile = P.C1; end
  end
  if nargin>3
    fprintf('Ignoring additional parameters\n');
  end
  SET = regexprep(maffile,'^.*/([^/])+$','$1');
else
  % standard parameter structure:   plotMutationSpectrumCategLegos(WORKSPACE,SET,MUTSIG,OUT,P)
  if (nargin<3)
    MUTSIG='mutsig1.5'
  end
  if (nargin<4)
    OUT=['~/' WORKSPACE '/plots/' datestr(now, 'ddmmmyyyy')];
  end
  
  slashjob='/';
  if ismember(P.FH_AREA,'/local/cga-fh/cga/')
      slashjob='/jobs/';
  end
  maffile=[P.FH_AREA '/' WORKSPACE '/Individual_Set/' SET slashjob  MUTSIG '/' SET '.final_analysis_set.maf'];
  covfile=[P.FH_AREA '/' WORKSPACE '/Individual_Set/' SET slashjob  MUTSIG '/' SET '.coverage.mat'];
  if isfield(P,'C1'), covfile = P.C1; end
end



if ~exist('OUT','var')
  OUT = '.';
end

if ~exist(OUT,'dir')
    fprintf('creating output area: %s\n',OUT);
    if (P.unix)
        unix(['mkdir -p ' OUT],'-echo');
    else
      error(' need output area ')
    end
end
if isstruct(P)
    if isfield(P,'FH_AREA') 
    	FH_AREA=P.FH_AREA;
    end
    if (isfield(P,'X') & isfield(P,'C1') )
        X=P.X;
        if isfield(X,'Variant_Type')
          vt = X.Variant_Type;
        elseif isfield(X,'classification')
          vt = X.classification;
        else
          error('help me!');
        end
        X=trimStruct(X,find(ismember(vt,'SNP')));
        C1=P.C1;
        fprintf('***\n override input with P.X and P.C1\n')
 	    fprintf(['MutationSpectrumCategs\n' ]);
        [XY,X,C1]=MutationSpectrumCategs(X,C1);

    elseif (isfield(P,'X') & isfield(P,'custom_WGS_covfile') )
        X=P.X;
        if isfield(X,'Variant_Type')
          vt = X.Variant_Type;
        elseif isfield(X,'classification')
          vt = X.classification;
        else
          error('help me!');
        end
        X=trimStruct(X,find(ismember(vt,'SNP')));
        C1=P.custom_WGS_covfile;
        fprintf('***\n override input with P.X and P.C1\n')
            fprintf(['MutationSpectrumCategs\n' ]);
        [XY,X,C1]=MutationSpectrumCategs(X,C1, [], true);
    else
      if (ischar(covfile))
          fprintf(['load\t' maffile '\nload\t' covfile '\nMutationSpectrumCategs\n' ]);
      else
          fprintf(['load\t' maffile '\nuser coverage \nMutationSpectrumCategs\n' ]);    
      end
      [XY,X,C1]=MutationSpectrumCategs(maffile,covfile);
    end

end

[q k]=sort(X.patient);
X=trimStruct(X,k);

fps=[OUT '/' SET '.mutation_legos.ps'];
fpdf=[OUT '/' SET '.mutation_legos.pdf'];
fhtml=[OUT '/' SET '.mutation_legos.html'];
%%


catfile='~/Projects/Cancer/tools/matlab/categs.txt';
if (~exist(catfile))
    %catfile='/xchip/cga1/lawrence/db/hg19/c65e/categs.txt';
    catfile='/xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/categs.txt';
end

categs=load_table(catfile);
k = map_categories_to_65(catfile); categs.longname=categs.name; q=regexp(categs.name,':','split');
categs.name=cellfun(@(x) x(1),q); [C.name i]=unique(categs.name);
C.id=k(i); C=trimStruct(C,sort(C.id)); C.N=length(C.id);

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
Np=length(unique(X.patient));
z=1e6;

if (round(Np*prod(size(XY.ncov)))==round(sum(XY.ncov(:))))
    z=1;
    zlab=' count per sample ';
    fprintf('mutation counts per sample \n');
else    
    fprintf('mutation rate lego\n');
    zlab='rate/Mb ';
end
%if (~P.visible), figure('Visible','Off'); end
clf
fp=get(gcf,'position');
set(gcf,'position',[fp(1:2) 924 604])

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
if (P.zscale), xp=0.15; end
text(xp,0.95,[SET ' '],'units','normalized','fontsize',14,'interpreter','none');
text(xp,0.9,sprintf('n=%d',X.N),'units','normalized','fontsize',14)    

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
	axes('Position', [0.59,0.07, .45, .5]);
else
	axes('Position', [0.59,0.05, .45, .5]);
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
f=[OUT '/' SET '.mutation_profile'];
f1=[SET '.mutation_profile.png'];
%print(gcf,'-depsc','-r400',[f '.eps'])
saveas(gcf,[f '.fig'],'fig')
saveas(gcf,[f '.png'],'png')
print(gcf,'-dill','-painters',[f '.ai'])
print(gcf,'-dpsc',fps);
fhtml=fopen(fhtml,'wt');
fprintf(fhtml,'<BODY>\n');
if exist('WORKSPACE','var') && exist('SET','var') && exist('MUTSIG','var')
  fprintf(fhtml,'<P> Workspace:\t %s \n',WORKSPACE);
  fprintf(fhtml,'<br> Individual_Set:\t %s \n',SET);
  fprintf(fhtml,'<br> MutSig version:\t %s\n',MUTSIG);
else
  fprintf(fhtml,'<br> MAF file:\t %s\n',maffile);
end
fprintf(fhtml,'<br>Mutation rate spectrum\n');
fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);


%%

fprintf('mutation count lego\n');
clf
bar3_with_colors(XY.n,XY.col)
set(gca,'xticklabel',[],'yticklabel',[],'ylim',[0.5 8.5]) %,'ztick',[],'visible','off');
set(gcf,'color',[1 1 1])
zlabel('counts ');

h0=gca;
text(0.05,0.95,[SET ' '],'units','normalized','fontsize',14,'interpreter','none');
text(0.05,0.9,sprintf('n=%d',X.N),'units','normalized','fontsize',14)    
q=regexp(X.CT,'-','split');
q=cellfun(@(x) x(2),q);
SNVT={'C>T';'C>A';'C>G';'A>G';'A>C';'A>T'};
nc=zeros(6,1);
for i=1:length(SNVT)
  k=strfindk(q(:),SNVT{i});
  nc(i)=length(k)+1e-5;
end
axes('Position', [0.05,0.6, .2, .2]);
h=pie(nc,repmat({''},length(nc),1));
colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7;];
for i=1:6
  set(h(2*i-1),'facecolor',colors(i,:))
end;
h1=legend(regexprep(SNVT,'>','>'));
set(h1,'Position',[.25 .6 .1 .2]);
set(h,'visible','off')
q=regexprep(XY.cat(1:4,1:4),'-','_'); 
box1=cellfun(@(x) [x([1 2 7]) ' '],q,'UniformOutput',false);
axes('Position', [0.59,0.05, .45, .5]);
h=bar3(ones(4,4),'w');
set(h,'EdgeAlpha',0.25);
set(gca,'zlim',[0 50]);
%bar3([0 1],[0 1],[0 1],'color','w')
for i=1:4,for j=1:4
  text((j-0.35),(i)+0.15,4,box1(i,j),'horizontalAlign','left','EraseMode','none','interpreter','none') %,'units','normalized')
end,end
set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);
set(h1,'OuterPosition',[0.7871 0.5546 0.1616 0.2670]);
set(h1,'Position', [0.7942 0.5594 0.1509 0.2622]);

%%
f=[OUT '/' SET '.mutation_profile_counts'];
%saveas(gcf,[f '.eps'],'epsc');
%print(gcf,'-depsc','-r400',[f '.eps'])
saveas(gcf,[f '.fig'],'fig');
print(gcf,'-dpsc','-append',fps);
saveas(gcf,[f '.png'],'png');
print(gcf,'-dill','-painters',[f '.ai'])
f1=[SET '.mutation_profile_counts.png'];
fprintf(fhtml,'<hr><P>Mutation count spectrum\n');
fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);



%%
clf
%figure
w=0; p0=1;
pat=unique(X.patient);
if isnumeric(pat)
   X.patient=cellstr(num2str(X.patient));
   pat=unique(X.patient);
end   
Npat=length(pat);
NX=4; NY=3; DX=1/(NX+0.7); DY=1/(NY+0.3);
    


if ~isfield(X,'i_tumor_f')
    if isfield(X,'t_alt_count')&isfield(X,'t_ref_count')
       X.i_tumor_f=X.t_alt_count./(X.t_alt_count+X.t_ref_count);
    else
       X.i_tumor_f=NaN(size(X.Hugo_Symbol));
    end
end
X.pat=zeros(size(X.i_tumor_f));
Z=[];
for p=1:Npat
    w=w+1;
    fprintf('%d of %d %s\n',p,Npat,pat{p})
    if (w>(NX*NY))
        f=[OUT '/' SET '.mutation_profile_samples_' num2str(p0) '-' num2str(p-1) ];
        saveas(gcf,[f '.fig'],'fig');
        %saveas(gcf,[f '.eps'],'epsc');
        %print(gcf,'-depsc','-r400',[f '.eps'])
        saveas(gcf,[f '.png'],'png');
        print(gcf,'-dpsc','-append',fps);
        print(gcf,'-dill','-painters',[f '.ai'])
        f1=[SET '.mutation_profile_samples_' num2str(p0) '-' num2str(p-1) '.png'];
        q1=1;
        fprintf(fhtml,'<br><hr>');
        fprintf(fhtml,'<table  width="80%%">\n<tr>\n');
        for p1=p0:(p-1)
            fprintf(fhtml,'<td>%d. %s</td>\n',p1,pat{p1});
            if (q1==4)
                fprintf(fhtml,'</tr><tr>');
                q1=0;      
            end
            q1=q1+1;    
        end
        fprintf(fhtml,'</tr>\n</table >\n');
        fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);
        w=1;
        clf;
        %figure;
        p0=p;
        
    end
    
    [ix,iy]=ind2sub([NX,NY],w); iy=NY-iy+1;
    
    zx=1+(P.zscale)*0.2;
        
    subplot('position',[0.07+(ix-1)*DX 0.1+(iy-1)*DY DX/zx DY])
   
    %subplot(5,2,w)
    %X1=trimStruct(X,strfindk(X.patient,pat{p}));
    k1=find(ismember(X.patient,pat{p}));
    X1=trimStruct(X,k1);
    C11=C1;
    C11.orig_cov=C11.orig_cov(p,:);
    %C11=trimStruct(C1,p)
    [XY1]=MutationSpectrumCategs(X1,C11);
    Q(p).name=pat{p};
    Q(p).n=XY1.n(:);
    Q(p).r=z*XY1.n(:)./(XY1.ncov(:));
    Q(p).c=XY1.cat(:);
    Q(p).NTOT=sum(XY1.n(:));
    Q(p).NART=sum(XY1.n(kART));
    Q(p).AF=median(X1.i_tumor_f);
    Q(p).FART=XY1.n(kART)/sum(XY1.n(:));
    Q(p).AF=median(X1.i_tumor_f);
    Q(p).AFART=NaN;
    Q(p).AFnART=NaN;
    Q(p).AFCpGT=NaN;
    k=find(ismember(X1.CT,{ART}));
    if (~isempty(k))
        Q(p).AFART=median(X1.i_tumor_f(k));
    end
    k=strfindk(X1.CT,ARTX,'v');
    if (~isempty(k))
        Q(p).AFnART=median(X1.i_tumor_f(k));
    end
    k=strfindk(X1.CT,CpGT);
    if (~isempty(k))
        Q(p).AFCpGT=median(X1.i_tumor_f(k));
    end
    X.pat(k1)=repmat(p,size(k1));
    bar3_with_colors(z*XY1.n./XY1.ncov,XY1.col)
    
    if (~P.zscale)
        set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);
    else
        set(gca,'xticklabel',[],'yticklabel',[],'ylim',[0.5 8.5],'fontsize',8);
    end

    %set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');
    set(gcf,'color',[1 1 1]);
    text(0.05,0.95,pat{p},'units','normalized','fontsize',round(14/sqrt(NY)),'interpreter','none')    
    text(0.05,0.85,sprintf('n=%d',X1.N),'units','normalized','fontsize',round(14/sqrt(NY)))  
    Z.name{p,1}=Q(p).name;
    Z.NTOT(p,1)=Q(p).NTOT;
    Z.NART(p,1)=Q(p).NART;
    Z.NnART(p,1)=sum(XY1.n(kARTX));
    Z.NCpGT(p,1)=length(strfindk(X1.CT,CpGT));
    Z.FART(p,1)=Q(p).FART;
    Z.FnART(p,1)=Z.NnART(p,1)/Z.NTOT(p,1);
    Z.FCpGT(p,1)=Z.NCpGT(p,1)/Z.NTOT(p,1);
    Z.AF(p,1)=Q(p).AF;
    Z.AFART(p,1)=Q(p).AFART;
    Z.AFnART(p,1)=Q(p).AFnART;
    Z.AFCpGT(p,1)=Q(p).AFCpGT;
end

f=[OUT '/' SET '.mutation_profile_samples_' num2str(p0) '-' num2str(p) ];
saveas(gcf,[f '.fig'],'fig');
%saveas(gcf,[f '.eps'],'epsc');
%print(gcf,'-depsc','-r400',[f '.eps'])
saveas(gcf,[f '.png'],'png');
print(gcf,'-dpsc','-append',fps)
print(gcf,'-dill','-painters',[f '.ai'])

f1=[SET '.mutation_profile_samples_' num2str(p0) '-' num2str(p) '.png'];
q1=1;
fprintf(fhtml,'<hr><table  width="80%%">\n<tr>\n');
for p1=p0:p
    fprintf(fhtml,'<td>%d. %s</td>\n',p1,pat{p1});
         
    if (q1==4)
        fprintf(fhtml,'</tr> ');
        fprintf(fhtml,'<tr> ');
        q1=0;      
    end
    q1=q1+1;    
end
fprintf(fhtml,'</tr>\n</table>\n');
fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);

%%
fprintf('allele frequency slice legos\n')
clf;
NX=2; NY=2; DX=1/(NX+0.1); DY=1/(NY+0.3);
%figure;
AF=[0 0.1 0.25 0.5 1];
for w=1:4
    [ix,iy]=ind2sub([NX,NY],w); iy=NY-iy+1;
    X1=trimStruct(X,find( (X.i_tumor_f>=AF(w))&(X.i_tumor_f<AF(w+1))) );
    if (X1.N<1), continue; end
    [XY1]=MutationSpectrumCategs(X1,C1);
    subplot('position',[0.05+(ix-1)*DX 0.1+(iy-1)*DY DX DY])   
    bar3_with_colors(z*XY1.n./XY1.ncov,XY1.col)
    s=sprintf('%.2f<=AF<%.2f',AF(w),AF(w+1));
    set(gca,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);
	if (w==1)
	    text(0.75,0.99,[SET ' '],'units','normalized','fontsize',14,'interpreter','none')    
	end
    text(0.05,0.87,s,'units','normalized','fontsize',round(14/sqrt(NY)),'interpreter','none')    
    text(0.05,0.80,sprintf('n=%d',X1.N),'units','normalized','fontsize',round(14/sqrt(NY)))    
end

f=[OUT '/' SET '.mutation_profile_AF' ];
saveas(gcf,[f '.fig'],'fig')
%saveas(gcf,[f '.eps'],'epsc')
%print(gcf,'-depsc','-r400',[f '.eps'])
saveas(gcf,[f '.png'],'png')
print(gcf,'-dpsc','-append',fps)
print(gcf,'-dill','-painters',[f '.ai'])

f1=[SET '.mutation_profile_AF.png'];
fprintf(fhtml,'<hr><P>Mutation count spectrum sliced by allele fraction\n');
fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);




%%
fprintf('allele frequency histogram figure\n')

%allele fraction histos  
clf
%ARTX=ART((2:end)-1);
k1=find(ismember(X.CT,{ART}));
k2=strfindk(X.CT,ARTX,'v');
k3=strfindk(X.CT,CpGT);
afb=0:0.01:1;
naf=0;
afart=NaN;
if (~isempty(k1))
    naf=hist(X.i_tumor_f(k1),afb);
    afart=median(X.i_tumor_f(k1));
end

%PY0=0.2;PX0=0.15; DX=0.8; DY=0.4
subplot(2,1,1)
%subplot('position',[PX0 PY0+0.55 DX DY])   
    
nx=hist(X.i_tumor_f(k2),afb);
afx=median(X.i_tumor_f(k2));
x1=[afb' afb']'; x1=x1(2:end);
y1=[naf' naf']'; y1=y1((2:end)-1);
y2=[nx' nx']'; y2=y2((2:end)-1);
h=plot(x1,y1,'k-',x1,y2,'b-',x1,0*x1,'k');
set(h,'linewidth',1.5);
if exist('WORKSPACE','var') && exist('SET','var')
  text(0.97,0.9,[WORKSPACE ' '],'units','normalized','interpreter','none','horizontalalign','right');
  text(0.97,0.85,[SET ' '],'units','normalized','interpreter','none','horizontalalign','right');
else
  text(0.97,0.9,[maffile ' '],'units','normalized','interpreter','none','horizontalalign','right');
end
text(0.97,0.8,sprintf('median=%.2f',afart),'units','normalized','interpreter','none','horizontalalign','right');
text(0.97,0.75,sprintf('median=%.2f',afx),'units','normalized','interpreter','none','horizontalalign','right','color','b');
hl=legend({'CCA>CAG','non-CC>CA  '},'location','E');
%text(0.95,0.7,sprintf('non-CC>CA mutations '),'color','b','units','normalized','interpreter','none','horizontalalign','right');
%text(0.95,0.6,sprintf('n=%d, median AF=%.2f ',sum(nx),afx),'color','b','units','normalized','interpreter','none','horizontalalign','right');
%text(0.95,0.5,sprintf('CCA>CAG mutations '),'color','k','units','normalized','interpreter','none','horizontalalign','right');
%text(0.95,0.4,sprintf('n=%d, median AF=%.2f ',sum(naf),afart),'color','k,'units','normalized','interpreter','none','horizontalalign','right');
if any(~isnan(X.i_tumor_f(k2)))
    axis([0 1 0 1.1*max([y1 y2])]);
end
textobj = findobj('type', 'text');
set([textobj; hl], 'fontunits', 'points','fontsize', 10);
xlabel('allele fraction ');
ylabel('mutations ');
set(gca,'xtick',0:0.1:1,'xticklabel',0:0.1:1)
%grid on


%tumor_lod histos  
subplot(2,1,2)
lodb=0:0.3:100;
nlod=0;
lodart=NaN;
olod=NaN;

if ~isfield(X,'i_t_lod_fstar')
    if ~isfield(X,'Start_position') && isfield(X,'start'), X = rename_field(X,'start','Start_position'); end
    X.i_t_lod_fstar=0*X.Start_position;
end
if (~isempty(k1))
    nlod=hist(X.i_t_lod_fstar(k1),lodb);
    olod=nlod(end); nlod(end)=0;
    lodart=median(X.i_t_lod_fstar(k1));
end

   
nx=hist(X.i_t_lod_fstar(k2),lodb);
ox=nx(end); nx(end)=0;
lodx=median(X.i_t_lod_fstar(k2));
x1=[lodb' lodb']'; x1=x1(2:end);
y1=[nlod' nlod']'; y1=y1((2:end)-1);
y2=[nx' nx']'; y2=y2((2:end)-1);
h=plot(x1,y1,'k-',x1,y2,'b-',x1,0*x1,'k');
set(h,'linewidth',1.5);
text(0.97,0.9,sprintf('median=%.1f, overflow=%d',lodart,olod),'units','normalized','interpreter','none','horizontalalign','right');
text(0.97,0.85,sprintf('median=%.1f, overflow=%d',lodx,ox),'units','normalized','interpreter','none','horizontalalign','right','color','b');

%hl=legend({'CCA>CAG','non-CC>CA  '},'location','E');
if any(X.i_t_lod_fstar(k2)>0)
    axis([min(lodb) max(lodb) 0 1.1*max([y1 y2])]);
end
textobj = findobj('type', 'text');
set([textobj; hl], 'fontunits', 'points','fontsize', 10);
xlabel('tumor_lod ');
ylabel('mutations ');
%set(gca,'xtick',0:5:50)
%grid on
f=[OUT '/' SET '.mutation_histogram' ];
saveas(gcf,[f '.fig'],'fig');
%saveas(gcf,[f '.eps'],'epsc');
%print(gcf,'-depsc','-r400',[f '.eps'])
saveas(gcf,[f '.png'],'png')
print(gcf,'-dpsc','-append',fps);
print(gcf,'-dill','-painters',[f '.ai'])
if (P.unix)
  [status, result]=unix(['ps2pdf ' fps ' ' fpdf]);
end
f1=[SET '.mutation_histogram.png'];
fprintf(fhtml,'<hr><P>Mutation allele fraction distribution \n');
fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);
fprintf(fhtml,'<br>Mutation tumor_lod \n');

clf
%subplot(2,1,1)
subplot('position',[0.1 0.4 0.9 0.5])   

q=regexp(X.CT,'-','split');
SNV=cellfun(@(x) x(2),q);
SNVT={'C>T';'C>A';'C>G';'A>G';'A>C';'A>T'};
NSNV=zeros(Npat,6);
NART=zeros(Npat,1);
NARTX=zeros(Npat,1);
for p=1:Npat
  kp=find(ismember(X.patient,pat{p}));
  SNV1=SNV(kp);
  [i1 k1]=ismember(SNV1,SNVT);
  for i=1:length(SNVT)
    NSNV(p,i)=sum(k1==i);
  end
  [i1 k1]=ismember(X.CT(kp),'C-C>A-G');  
  NART(p)=sum(i1);
  k1=strfindk(X.CT(kp),'C-C>A');  
  NARTX(p)=length(k1);
end

c0= 0.25*[1 1 1 ];
colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7;];
h=bar(sum(NSNV,2),'facecolor',colors(1,:),'edgecolor', c0);
set(h,'BaseValue',0.9)
set(gca,'yscale','log');
YMAX=max(sum(NSNV'));
% put CT bins at bottom
NCT=NSNV(:,2);
NSNV(:,2)=0;
hold on
for t=2:6
    h(t)=bar(sum(NSNV(:,t:end),2),'facecolor',colors(t,:),'edgecolor', c0);
end
bar(0*NARTX,'facecolor',c0,'edgecolor', c0);
bar(0*NART,'facecolor','k','edgecolor', c0);
h(7)=bar(NCT,'facecolor',colors(2,:),'edgecolor', c0);
h(8)=bar(NARTX,'facecolor',c0,'edgecolor', c0);
h(9)=bar(NART,'facecolor','k','edgecolor', c0);
hold off
axis([0.5 Npat+0.5 0.9 1.5*YMAX])
%set(gca,'ylim',[0 1.1*max(sum(NSNV,2))]);
textobj = findobj('type', 'text');
set(textobj, 'fontunits', 'points');
set(textobj, 'fontsize', 12);
ylabel('mutations ');
SL=[SNVT; 'CCX>CAX'; 'CCG>CAG' ];
SL(2)={'C>A,no CC'};
hl=legend(SL,'location','EO');
textobj = findobj('type', 'text');
set([textobj; hl], 'fontunits', 'points','fontsize', 10);
set(gca,'ytick',10.^(0:4))
if (Npat>40)
    xlabel('patients (ordered by name) ');
else
    set(gca,'xticklabel',[])
    xt=1:Npat;
    fsz=round(10*30/Npat); if (fsz>11), fsz=11; end
    text(xt,0.75+0*xt,pat,'rotation',-25,'fontsize',fsz)
end
grid on

%%
%pos0=get(gcf,'position')
%posW=pos0; posW(3:4)=[1500 500]; 
%set(gcf,'position',posW)
f=[OUT '/' SET '.mutation_CCG_histogram' ];
saveas(gcf,[f '.fig'],'fig');
%saveas(gcf,[f '.eps'],'epsc');
%print(gcf,'-depsc','-r400',[f '.eps'])
saveas(gcf,[f '.png'],'png')
print(gcf,'-dpsc','-append',fps);
print(gcf,'-dill','-painters',[f '.ai'])

if (P.unix)
  [status, result]=unix(['ps2pdf ' fps ' ' fpdf]);
end
f1=[SET '.mutation_CCG_histogram.png'];
fprintf(fhtml,'<hr><P>Mutation breakdown per sample \n');
fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);


%set(gcf,'position',pos0)

%%
fprintf('coverage\n');
clf
imagesc((XY.ncov)/ (Npat*1e6),[0 1.1*max(XY.ncov(:))]/(Npat*1e6)); b=colorbar; colormap(1-bone(100))
set(gca,'xtick',1:12,'xticklabel',repmat({'x_T','x_C','x_A','x_G'},1,3))
set(gca,'ytick',1:8,'yticklabel',{'TC_','CC_','AC_','GC_','TA_','CA_','AA_','GA_'})

xlabel('trailing base ');
ylabel('leading+center base ');
ylabel(b,'Mb bases covered per sample ');
axis equal
set(gca,'xlim',[0.5 12.5],'ylim',[0.5 8.5])
rectangle('position',[0.55 0.55 3.9 3.9],'edgecolor',XY.col(1,1,:),'linewidth',6)
rectangle('position',[0.55 4.55 3.9 3.9],'edgecolor',XY.col(5,1,:),'linewidth',6)
rectangle('position',[4.55 0.55 3.9 3.9],'edgecolor',XY.col(1,5,:),'linewidth',6)
rectangle('position',[4.55 4.55 3.9 3.9],'edgecolor',XY.col(5,5,:),'linewidth',6)
rectangle('position',[8.55 0.55 3.9 3.9],'edgecolor',XY.col(1,9,:),'linewidth',6)
rectangle('position',[8.55 4.55 3.9 3.9],'edgecolor',XY.col(5,9,:),'linewidth',6)
title('base coverage ')
text(1,1,'C>G','BackgroundColor',[1 1 1],'fontsize',14,'edgecolor','k')
text(1,5,'A>T','BackgroundColor',[1 1 1],'fontsize',14,'edgecolor','k')
text(5,1,'C>A','BackgroundColor',[1 1 1],'fontsize',14,'edgecolor','k')
text(5,5,'A>C','BackgroundColor',[1 1 1],'fontsize',14,'edgecolor','k')
text(9,1,'C>T','BackgroundColor',[1 1 1],'fontsize',14,'edgecolor','k')
text(9,5,'A>G','BackgroundColor',[1 1 1],'fontsize',14,'edgecolor','k')

f=[OUT '/' SET '.mutation_coverage'];
%saveas(gcf,[f '.eps'],'epsc');
%print(gcf,'-depsc','-r400',[f '.eps'])
saveas(gcf,[f '.fig'],'fig');
print(gcf,'-dpsc','-append',fps);
saveas(gcf,[f '.png'],'png');
f1=[SET '.mutation_coverage.png'];
fprintf(fhtml,'<hr><P>Mutation coverage spectrum\n');
fprintf(fhtml,'<img border="0" src="%s" type="image/png" alt="%s" width="924" height="604" />\n',f1,SET);

%%
fprintf(fhtml,'</BODY>\n');
fclose(fhtml);



%%

AF=median(X.i_tumor_f);
%artifact bin
AFART=NaN;
k=find(ismember(X.CT,{ART}));
if (~isempty(k))
    AFART=median(X.i_tumor_f(k));
end
%non-artifact bins
AFnART=NaN;
k=strfindk(X.CT,ARTX,'v');
if (~isempty(k))
    AFnART=median(X.i_tumor_f(k));
end
AFCpGT=NaN;
k=strfindk(X.CT,CpGT);
if (~isempty(k))
   AFCpGT=median(X.i_tumor_f(k));
end


fid=1;

f=[OUT '/' SET '.mutation_spectrum.txt'];
fid=fopen(f,'wt');
Nc=length(XY.n(:));
fmt1='%s\t%s\t%d\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f';
fmt2='%s\t%s\t%d\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f';
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s','name','#/Mb','Nindiv','OK','N_SNV','N_CCG>CAG','F_CCG>CAG','AF','AF_CCG>CAG','AF_nCC>CA','AF_CpG>T');
for i=1:Nc
    fmt1=[fmt1 '\t%d'];
    fmt2=[fmt2 '\t%.2f'];
    fprintf(fid,'\t%s',XY.cat{i});
end
fmt1=[fmt1 '\n'];
fmt2=[fmt2 '\n'];
fprintf(fid,'\n');

OK='OK';
if (fART>=0.2), OK='FAIL';end
if ((fART>0.05)&(fART<0.2)), OK='MAYBE';end

fprintf(fid,fmt1,SET,'#',Npat,OK,X.N,nART,fART,AF,AFART,AFnART,AFCpGT,XY.n(:));
if (P.printrates)
    fprintf(fid,fmt2,SET,'mutations/Mb',Npat,OK,X.N,nART,fART,AF,AFART,AFnART,AFCpGT,z*XY.n(:)./XY.ncov(:));
end

for p=1:Npat
    OK='OK';
    if (Q(p).FART>0.2), OK='FAIL';end
    if ((Q(p).FART>0.05)&(Q(p).FART<0.2)), OK='MAYBE';end

    fprintf(fid,fmt1,Q(p).name,'#',1,OK,Q(p).NTOT,Q(p).NART,Q(p).FART,Q(p).AF,Q(p).AFART,Q(p).AFnART,Q(p).AFCpGT,Q(p).n(:));
    if (P.printrates)
        fprintf(fid,fmt2,Q(p).name,'mutations/Mb',1,OK,Q(p).NTOT,Q(p).NART,Q(p).FART,Q(p).AF,Q(p).AFART,Q(p).AFnART,Q(p).AFCpGT,Q(p).r(:));
    end
end


fclose(fid);
fprintf('\n\n****\ndone \n****\n')

function progressreport(n,N,msg)
fprintf(1,'\b%d\t%d\t%s',n,N,msg); pause(.1)


function test

clear
cd ~/Projects/Cancer/tools/matlab
addpath('~/CancerGenomeAnalysis/trunk/matlab')
addpath('~/CancerGenomeAnalysis/trunk/matlab/seq')
addpath('~/CancerGenomeAnalysis/trunk/matlab/mike')


close all
%WORKSPACE='TCGA_Bladder' ; SET='BLCA53-new'
WORKSPACE='Sigma_Pediatric';MUTSIG='mutsig_1.5_WES'
SET='AN_SIGMA_Rhabdoid_16Feb2012_Diagnostic'
%SET='AN_SIGMA_Rhabdoid_16Feb2012_TX'
%SET='AN_SIGMA_EwingsSarcoma_20Dec2011'
%SET='PR_SIGMA_Rhabdoid_Capture'
%WORKSPACE='CLL_paper_no_cutoff'; SET='An_CLL_February2012_v4';MUTSIG='mutsig1.5'
%WORKSPACE='An_Neuroblastoma_Paper_Analysis'; 
%SET='An_Neuroblastoma_Genome_Exome_99'; MUTSIG='mutsig_1.5'
%WORKSPACE='AN_THCA'
%SET='AN_TCGA_THCA_22Feb2012_AB'
%MUTSIG='mutsig1.5'

%TODAY=datestr(now, 'ddmmmyyyy')
%OUT=['~/' WORKSPACE '/plots/' TODAY];
[XY,X,C1]=plotMutationSpectrumCategLegos(WORKSPACE,SET,MUTSIG)

function test_eso_wgs
clear
cd ~/Projects/Cancer/ESO
WORKSPACE='An_ESO'
SET='PR_Esophageal_CIP_WGS'
MUTSIG='mutsig1.5'
OUT=['~/Projects/lego_plots/'  WORKSPACE '/' datestr(now, 'ddmmmyyyy') '_AVECOV'];
unix(['mkdir -p ' OUT])
P.C1='genome'
P.zscale = true;
P.printrates=true;
set(gcf,'renderer','painters');
f='/local/cga-fh/cga/Test_OxoG/Individual_Set/PR_Esophageal_CIP_WGS/jobs/wgs/annotation_workflow/aggregated_snp_maf/PR_Esophageal_CIP_WGS.snp.maf.annotated'
unix(['cut -f1,5-17,64-65,81-86 ' f ' > PR_Esophageal_CIP_WGS.snp.maf'])
f='PR_Esophageal_CIP_WGS.snp.maf'
P.X=load_table(f)
[XY,X,C1]=plotMutationSpectrumCategLegos(WORKSPACE,SET,MUTSIG,OUT,P)

q=load('eso_WGS_coverage.mat')
Npat=length(unique(P.X.Tumor_Sample_Barcode));
P.C1.orig_cov=repmat(q.N',Npat,1);
OUT=['~/Projects/lego_plots/'  WORKSPACE '/' datestr(now, 'ddmmmyyyy') '_ESOCOV'];
[XY,X,C1]=plotMutationSpectrumCategLegos(WORKSPACE,SET,MUTSIG,OUT,P)

