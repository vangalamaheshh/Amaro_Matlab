clear
% capSeq Seg file
DATA=1% 'capseg'
DATA=2% 'segseq'
DATA=3% 'segseq primaries'
DATA=4% NSCLC
DATA=5% PrimaryOnly27SampleSegments_Gistic
if DATA==1
    AC='/local/cga-fh/cga/An_MOUSE_SCLC_PAIR/Pair_Set/'
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/Mouse_OriginalSet_Primary/'
    SET='Original_TP'
    A1=[AC SET '/jobs/CapSegConcatenateSegFiles/']
    f1='seg_file.seg_file'
elseif DATA==2
    AC='/local/cga-fh/cga/An_MOUSE_SCLC_PAIR/Pair_Set/'
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/WGS_Gistic/'
    SET='LowPass'
    A1='/Users/stewart/Projects/Cancer/mouse/SegSeq/LP/'
    f1='lp.18Feb2013.seg.txt'
elseif DATA==3
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/WGS_Gistic/'
    SET='LowPass.Primary'
    A1='/Users/stewart/Projects/Cancer/mouse/Gistic2/LowPassPrimaries_only/Primaries_only/'
    f1='PrimaryWGS_lowpass.seg'
elseif DATA==4   
    AC='/Users/stewart/Projects/Cancer/mouse/Gistic2/NSCLC/'
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/NSCLC/'
    SET='NSCLC'
    A1='/Users/stewart/Projects/Cancer/mouse/Gistic2/NSCLC/'
    f1='NSCLC.seg.txt'    
elseif DATA==5
%     AC='/Users/stewart/Projects/Cancer/mouse/Gistic2/PrimaryOnly27SampleSegments_Gistic/'
%     AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/PrimaryOnly27SampleSegments_Gistic/'
%     SET='PrimaryOnly27SampleSegments_Gistic'
%     A1='/Users/stewart/Projects/Cancer/mouse/Gistic2/PrimaryOnly27SampleSegments_Gistic/'
%     f1='PrimaryOnly27SampleSegments_Gistic.seg.txt' 
%     
     AC='/Users/amaro/Downloads/PrimaryOnly27SampleSegments_Gistic/'
    AG='/Users/amaro/Downloads/PrimaryOnly27SampleSegments_Gistic/'
    SET='PrimaryOnly27SampleSegments_Gistic'
    A1='/Users/amaro/Downloads/PrimaryOnly27SampleSegments_Gistic/'
    f1='PrimaryOnly27Samples.seg' 
    %f1='NSCLC.seg.txt'    
end


VS=load_table([A1 f1]);
f2='broad_values_by_arm.txt'
VA=load_table([AG f2]);
f4='focal_data_by_genes.txt'
QF=load_table([AG f4]);
f5='broad_data_by_genes.txt'
QB=load_table([AG f5]);

XSCN=[];
ff=fieldnames(QF); s=ff(strfindk(ff,'SCLC'))'; 
if DATA==4   
     s=ff(strfindk(ff,'DM'))'; 
     VS.Sample=regexprep(VS.Sample,'\.','_')
end

XSCN.ff=s;
XSCN.sample=regexprep(s,'_seg','');
XSCN.gene=QF.Gene_Symbol;
XSCN.Cytoband=QF.Cytoband;
XSCN.ARM=cellfun(@(x) x{1},regexp(XSCN.Cytoband,'\w+[p|q]','match'),'UniformOutput',false);
ARMS=unique(sort(XSCN.ARM));
NS=length(XSCN.sample);
NG=length(XSCN.gene);
ff=fieldnames(QF);
XSCN.focal=NaN(NG,NS);
XSCN.broad=NaN(NG,NS);
XSCN.arm=NaN(NG,NS);
if (DATA==1)
    XSCN.sample=regexprep(XSCN.sample,'_','-');
end
for i=1:NS
    if (DATA==1)
        s=regexprep(XSCN.sample(i),'-','_');
        s=regexprep(s,'_Tumor','');    
    else
        s=XSCN.ff;
    end
    if ismember(s,ff)
        XSCN.focal(:,i)=QF.(s{1});
        XSCN.broad(:,i)=QB.(s{1});
        for a=1:length(ARMS)
            k=find(ismember(XSCN.ARM,ARMS{a}));
            XSCN.arm(k,i)=median(XSCN.broad(k,i));
        end
    end
end
clear QF QB

save(['./' SET '.SCNA.mat'],'-struct','XSCN')
save(['./' SET '.SCNA.structs.mat'],'XSCN','VS','VA')

%%

MM9=load_table('/Users/amaro/Downloads/MouseCytobandFile.txt')
Nc=length(MM9.cytoband)
MM9.arm=cellstr(cellfun(@(x) x{1},regexp(MM9.cytoband,'[p-q]','match'),'UniformOutput', false))
MM9.arm=strcat(regexprep(MM9.chrn,'chr',''),'q')
MM9A=[];
a=unique(MM9.arm);
for i=1:length(a)
    k=find(ismember(MM9.arm,a{i}));
    MM9A.arm(i,1)=a(i);
    MM9A.chrn(i,1)=MM9.chrn(k(1));
    MM9A.start(i,1)=min(MM9.start(k));
    MM9A.end(i,1)=max(MM9.end(k));
end
MM9A.x1=xhg19(chrom2num(MM9A.chrn),MM9A.start,'mm9');
MM9A.x2=xhg19(chrom2num(MM9A.chrn),MM9A.end,'mm9');
[q k]=sort(MM9A.x1);
MM9A=trimStruct(MM9A,k)
printStruct(MM9A,[],'arms.txt')


%%

clear
DATA=4  % NSCLC
DATA=5  % PrimaryOnly27SampleSegments_Gistic

if (DATA==1)
    SET='Original_TP'
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/Mouse_OriginalSet_Primary/'
elseif DATA==2
    SET='LowPass'
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/WGS_Gistic/'
elseif DATA==3
    SET='LowPass.Primary'
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/WGS_Gistic/'
elseif DATA==4
    SET='NSCLC'
    AG='/Users/stewart/Projects/Cancer/mouse/Gistic2/NSCLC/'
elseif DATA==5
    AG='/Users/amaro/Downloads/PrimaryOnly27SampleSegments_Gistic/'
    SET='PrimaryOnly27SampleSegments_Gistic'    
end
load(['./' SET '.SCNA.structs.mat'])
XSCN.sample=regexprep(XSCN.sample,'MOUSE-S','MOUSE_S')
f=[AG 'scores.gistic']
G=load_table(f)
G.x1=xhg19(G.Chromosome,G.Start,'mm9')/1e6
G.x2=xhg19(G.Chromosome,G.End,'mm9')/1e6
G.nlog10q=G.v__log10_q_value_;
Q=load([AG 'all_thresholded.by_genes.mat'],'T')
Q.B=load([AG 'broad_results.mat'])
f=['arms.txt']
MM9A=load_table(f)
f='/Users/amaro/Downloads/PR_MOUSE_SCLC_PAIR_Capture_All_Pairs-1.purity.11Dec2012.tsv'
P=load_table(f)
P.sample_id=cellstr(cellfun(@(x) [x{1} '-' x{2} '-Tumor-' x{5} '-' x{6}],regexp(P.pair,'-','split'),'UniformOutput', false))
%addpath('./CancerGenomeAnalysis/trunk/matlab/snp')
TX=load_table('/Users/amaro/Downloads/mouse_mutation_summary_dgm.txt')
if (DATA==1)
    Z=load_table('~/Projects/Cancer/mouse/FH/An_MOUSE_SCLC_PAIR.Original_TP.samples.txt')
    Z=trimStruct(Z,find(cellfun(@length,strfind(Z.sample_id,'-Tumor'))))
    Z.sample=regexprep(Z.sample_id,'-Tumor','')
    Z.individual_id=regexprep(Z.individual_id,'_','-')

elseif (DATA==5)
    Z=load_table('/Users/amaro/Documents/sample_ids.txt')
    Z=trimStruct(Z,find(cellfun(@length,strfind(Z.sample_id,'-Tumor'))))
    Z.sample=regexprep(Z.sample_id,'-Tumor','')
    Z.individual_id=regexprep(Z.individual_id,'_','-')
    s=fieldnames(VA); s=s(4:end);s=cellfun(@(x) [x{1} '_' x{2} '-' x{3} '-' x{6} '-' x{7}],regexp(s,'_','split'),'UniformOutput',false )
    [i m]=ismember(Z.sample,s)
    Z=trimStruct(Z,find(i))
elseif (DATA>4)
    Z=load_table('~/Projects/Cancer/mouse/FH/An_MOUSE_SCLC_PAIR.LowPass.samples.txt')
    Z=trimStruct(Z,strfindk(Z.sample_id,'-Normal','v'))
    Z=trimStruct(Z,strfindk(Z.sample_id,'-Tail','v'))
    Z.sample=regexprep(Z.sample_id,'MOUSE_SCLC_','MOUSE_SCLC-')
    Z.individual_id=cellstr(cellfun(@(x) [x{1} '-' x{2}  ],regexp(Z.sample,'-','split'),'UniformOutput', false))
    TX1.Tumor={ 'MOUSE_SCLC-1867-C1-TM'
        'MOUSE_SCLC-1867-T2-TP'
        'MOUSE_SCLC-1867-T3-TP'
        'MOUSE_SCLC-1868-N1-TM'
        'MOUSE_SCLC-3384-T1-TP'
        'MOUSE_SCLC-3384-T2-TP'
        'MOUSE_SCLC-3882-T5-TP'}
    TX1.Tumor=[TX1.Tumor; TX.Tumor]
    TX1.Bio_ID={ '1867C1'
        '1867T2'
        '1867T3'
        '1868N1'
        '3384T1'
        '3384T2'
        '3882T5'}
    TX1.Bio_ID=[TX1.Bio_ID; TX.Bio_ID]
    TX=TX1;
end
if (DATA<4)
    Z.Sample_label=Z.individual_id;
    [i m]=ismember(Z.sample,TX.Tumor)
    Z.sample_id(~i)
    Z.Sample_label(i)=TX.Bio_ID(m(m>0));
    [Z.individual_id Z.Sample_label]
    Z.Sample_label=regexprep(Z.Sample_label,'^AD','')
    Z.Sample_label=regexprep(Z.Sample_label,'^TP','')
    XSCN.sample_label=repmat({''},size(XSCN.sample))
    if DATA==1
        XSCN.sample_id=cellstr(cellfun(@(x) [x{1} '-' x{2} '-Tumor-' x{5} '-' x{6}],regexp(XSCN.sample,'-','split'),'UniformOutput', false))
        [i m]=ismember(XSCN.sample_id,Z.sample_id)
        XSCN.sample_label(i)=Z.Sample_label(m(m>0))
        XSCN.individual_id(i) = regexprep(Z.individual_id(m(m>0)),'MOUSE-SCLC-','')
    else
        XSCN.sample_label=cellfun(@(x) x{2},regexp(XSCN.sample,'_','split'),'UniformOutput', false)
        [i m]=ismember(XSCN.sample_label,Z.Sample_label)
        XSCN.sample_id=Z.sample_id(m(m>0))'
        q=char(XSCN.sample_label'); q1=q(:,1:4); q1=cellstr(q1)
        XSCN.individual_id=q1';
        q2=q(:,5:6); q2=cellstr(q2)
        
    end
elseif (DATA==5)
    XSCN.sample_id=regexprep(XSCN.sample,'_','-');
    XSCN.individual_id=cellfun(@(x) x(3),regexp(XSCN.sample_id,'-','split'))
    XSCN.sample_label=cellstr(cellfun(@(x) [x{3} '-' x{7}],regexp(XSCN.sample_id,'-','split'),'UniformOutput',false))

    
else
    XSCN.sample_label=regexprep(XSCN.sample,'_','-');
    XSCN.sample_id=regexprep(XSCN.sample,'_','-');
    q=char(XSCN.sample); q=q(:,1:6); 
    XSCN.individual_id=cellstr(q);
    
end

[XSCN.sample_id'  XSCN.sample_label' XSCN.individual_id']

[q k]=sort( XSCN.sample_label)

XSCN.sample=XSCN.sample(k)
XSCN.sample_id=XSCN.sample_id(k)
XSCN.individual_id=XSCN.individual_id(k)
XSCN.sample_label=XSCN.sample_label(k)
XSCN.focal=XSCN.focal(:,k)
XSCN.broad=XSCN.broad(:,k)
XSCN.arm=XSCN.arm(:,k)

%%

clf

if ~isfield(VS,'Segment_Mean')
    VS.Segment_Mean=VS.Copy_ratio;
end
VS1=trimStruct(VS,find(abs(VS.Segment_Mean)>0.2))

VS1=trimStruct(VS1,find(~ismember(VS1.Chromosome,'Y')));
subplot('position',[0.1 0.15 0.5 0.8])
if (~isnumeric(VS1.Chromosome))
    VS1.x1=xhg19(chrom2num(VS1.Chromosome),VS1.Start,'mm9')
    VS1.x2=xhg19(chrom2num(VS1.Chromosome),VS1.End,'mm9')
else
    VS1.x1=xhg19(VS1.Chromosome,VS1.Start,'mm9')
    VS1.x2=xhg19(VS1.Chromosome,VS1.End,'mm9')    
end
if (DATA==2)
    VS1.Sample=regexprep(VS1.Sample,'\.','-');
    VS1.Sample_id=regexprep(VS1.Sample,'-LP_ALL','')
    VS1.Sample_id=regexprep(VS1.Sample_id,'-NT$','')
    VS1.Sample_id=regexprep(VS1.Sample_id,'-TP-NT-','-Tumor-')
    VS1.Sample_id=regexprep(VS1.Sample_id,'-TP-X-','-Tumor-')
    VS1.Sample_id=regexprep(VS1.Sample_id,'-X$','')
    unique(sort(VS1.Sample_id))
    sort(XSCN.sample_id')
    [i m]=ismember(VS1.Sample_id,XSCN.sample_id);
    VS1=trimStruct(VS1,find(i))
elseif (DATA==3)
    unique(sort(VS1.Sample))
    XSCN.sample
    VS1.Sample_id=regexprep(VS1.Sample,'-','_')
    VS1.Sample_id=regexprep(VS1.Sample_id,'_seg','')
    unique(sort(VS1.Sample_id))
    XSCN.sample
    [i m]=ismember(VS1.Sample_id,XSCN.sample);
    VS1=trimStruct(VS1,find(i))     
elseif (DATA==4)
    unique(sort(VS1.Sample))
    XSCN.sample
    [i m]=ismember(VS1.Sample,XSCN.sample);
    VS1=trimStruct(VS1,find(i))  
elseif (DATA==5)
    VS1.Sample=regexprep(VS1.Sample,'\.','_')
    unique(sort(VS1.Sample))
    XSCN.sample
    [i m]=ismember(VS1.Sample,XSCN.sample);
    VS1=trimStruct(VS1,find(i))      
end


%XSCN.sample=regexprep(XSCN.sample,'-','_');
NS=length(XSCN.sample);
axis([-1 NS+1 0 max(VS1.x2)/1e6]); %set(gca,'visible','off')
set(gca,'box','on');
dx=0.99/2;
cg=gray(50);
c1=cg; c1(:,3)=1;
c2=flipud(cg); c2(:,1)=1;
c=[c1;c2];
N=length(VS1.x1)


for s=1:NS
    if (DATA<3)
        k=find(ismember(VS1.Sample_id,XSCN.sample_id{s}));
    elseif (DATA==3)
        k=find(ismember(VS1.Sample_id,XSCN.sample{s}));
    elseif (DATA>=4)
        k=find(ismember(VS1.Sample,XSCN.sample{s}));
    end
    x1=repmat([s-dx s+dx s+dx s-dx s-dx],length(k),1);
    yL=VS1.x1(k)/1e6;
    yH=VS1.x2(k)/1e6;
    y1=[yL yL yH yH yL];
    z1=VS1.Segment_Mean(k);
    if isfield(XSCN,'purity')
        SCN=(2*2.^VS1.Segment_Mean(k))
        SCN=2+((SCN-2)/XSCN.purity(s))
        SCN(SCN<0)=0.1;
        z1=log2(SCN/2);
    end
    kC=round(50*(z1+1)); kC(kC<1)=1; kC(kC>100)=100;
    c1=zeros(5,length(kC),3);
    for j=1:3
        for i=1:5
            c1(i,1:length(kC),j)=c(kC,j);
        end
    end
    patch(x1',y1',c1,'edgecolor','none')   
end


xe=xhg19(1:24,0*(1:24),'mm9')/1e6
xm=(xe(2:end)+xe((2:end)-1))/2
xe(end)=[];
NC=23
line([0;(NS+1)]*ones(1,NC),[1;1]*xe,'color',0.5*[1 1 1],'linestyle','--')
set(gca,'ydir','reverse')
%xlim([-1 NS+1.1])
xlim([-0 NS+0.5])
dx=0.5
xe=xhg19(1:24,0*(1:24),'mm9')/1e6
for i=[1:2:19]
    s=-0.0
    x1=[s-dx s+dx s+dx s-dx s-dx];
    y1=[xe(i) xe(i) xe(i+1) xe(i+1) xe(i)]'
    patch(x1,y1,'k','edgecolor','none')   
    s=NS+1
    x1=[s-dx s+dx s+dx s-dx s-dx];
    %patch(x1,y1,'k','edgecolor','none')      
end
set(gca,'xticklabel',[],'xtick',1:NS)
k=[1:12 13:2:19 23]
cL=num2chrom(k)
set(gca,'ytick',xm(k),'yticklabel',cL)
ylabel('Chromosome ')


uid=unique(sort(XSCN.individual_id))
mid=uid(1:2:end)
for m=1:2:length(uid)
    k=find(ismember(XSCN.individual_id,uid{m}))
    line(min(k)-0.5*[1 1],[-150 xe(end)],'linestyle','--','color',0.5*[1 1 1])
    line(max(k)+0.5*[1 1],[-150 xe(end)],'linestyle','--','color',0.5*[1 1 1])
    line([min(k)-0.5 max(k)+0.5],[10 10],'linestyle','-','linewidth',3,'color',0.5*[1 1 1])
end
set(gca,'xtick',1:NS,'xticklabel',XSCN.sample_label)
hl=rotateticklabel(gca,-40,-13.5)
line([-0.5 NS+0.5],[1 1],'linestyle','-','linewidth',1,'color',0*[1 1 1])
line([-0.5 NS+0.5],[1 1]*xe(end)-1,'linestyle','-','linewidth',1,'color',0*[1 1 1])
set(hl,'FontSize',12)

yL=[1 1.5 2 3 4]
yT=100*(log2(yL/2)+1)/2; yT(1)=5;
yL=regexp((num2str(yL,'%g\n')),'\n','split')
h=gca;
ylim([0 xe(end)])
%
axes('position',[.1  0.96  .2  .01])
image(reshape(c,1,100,3))
hi=gca;
set(hi,'ytick',[])
set(hi,'xtick',yT,'xticklabel',yL,'XAxisLocation','top')
text(-17,1,'SCN','fontsize',12)

xe=xhg19(1:24,0*(1:24),'mm9')/1e6
NC=length(xe)

%%
% broad
subplot('position',[0.62 0.15 0.08 0.8])
[i m]=ismember(Q.B.names,MM9A.arm)
xd=[-log10(Q.B.qD(i)+2e-12) -log10(Q.B.qD(i)+2e-12)]'; xd=xd(:);
yd=[MM9A.x1(m(m>0)) MM9A.x2(m(m>0))]'; yd=yd(:)/1e6;
xd(end+1)=0;
yd(end+1)=yd(end);
xd(xd>10)=10;
stairs(xd,yd,'b')
set(gca,'xdir','reverse','ydir','reverse')
xlim([-1 11])
k=[1:12 13:2:19 23]
set(gca,'ytick',xm(k),'yticklabel',[])
dx=0.95; NC=length(xe)
for i=1:2:19
    s=-1
    x1=[s-dx s+dx s+dx s-dx s-dx]';
    y1=[xe(i) xe(i) xe(i+1) xe(i+1) xe(i)]'
    patch(x1,y1,'k','edgecolor','none')   
end
ylim([0 xe(end)])
line([0;11]*ones(1,NC),[1;1]*xe,'color',0.5*[1 1 1],'linestyle','--')
hold on;
stairs(xd,yd,'b','linewidth',2)
hold off;
xt=get(gca,'xtick');
xt(1)=[];
set(gca,'xtick',xt);
line([0 11],[1 1],'linestyle','-','linewidth',1,'color',0*[1 1 1])
line([0 11 ],[1 1]*xe(end)-1,'linestyle','-','linewidth',1,'color',0*[1 1 1])
text(5,50,'del')

subplot('position',[0.7 0.15 0.08 0.8])
[i m]=ismember(Q.B.names,MM9A.arm)
xa=[-log10(Q.B.qA(i)+2e-12) -log10(Q.B.qA(i)+2e-12)]'; xa=xa(:);
ya=[MM9A.x1(m(m>0)) MM9A.x2(m(m>0))]'; ya=ya(:)/1e6;
xa(xa>10)=10;

stairs(xa,ya,'r')
xlim([-.05 11])
set(gca,'ytick',xm(k),'yticklabel',[])
ylim([0 xe(end)])
set(gca,'ydir','reverse')
line([0;11]*ones(1,NC),[1;1]*xe,'color',0.5*[1 1 1],'linestyle','--')
hold on;
stairs(xa,ya,'r','linewidth',2)
hold off;
text(0,-70,'Arm-level','fontsize',12,'horizontalAlign','center')
text(0,xe(end)+150,'-log10(q-value)','fontsize',12,'horizontalAlign','center')
line([0 11],[1 1],'linestyle','-','linewidth',1,'color',0*[1 1 1])
line([0 11 ],[1 1]*xe(end)-1,'linestyle','-','linewidth',1,'color',0*[1 1 1])
text(2.5,50,'gain')


% focal 
subplot('position',[0.80 0.15 0.08 0.8])
kd=find(ismember(G.Type,'Del'));
xd=[G.nlog10q(kd) G.nlog10q(kd)]'; xd=xd(:);
yd=[G.x1(kd) G.x2(kd)]'; yd=yd(:);
plot(xd,yd,'b','linewidth',2)
set(gca,'xdir','reverse','ydir','reverse')
xlim([-1 11])
set(gca,'ytick',xm(k),'yticklabel',[])

%xe=xhg19(1:24,0*(1:24),'hg19')/1e6
dx=0.9; 
for i=1:2:19
    s=-1
    x1=[s-dx s+dx s+dx s-dx s-dx]';
    y1=[xe(i) xe(i) xe(i+1) xe(i+1) xe(i)]'
    patch(x1,y1,'k','edgecolor','none')   
end
ylim([0 xe(end)])
line([0;11]*ones(1,NC),[1;1]*xe,'color',0.5*[1 1 1],'linestyle','--')
xt=get(gca,'xtick');
xt(1)=[];
set(gca,'xtick',xt);
ylim([0 xe(end)])
line([0 11],[1 1],'linestyle','-','linewidth',1,'color',0*[1 1 1])
line([0 11 ],[1 1]*xe(end)-1,'linestyle','-','linewidth',1,'color',0*[1 1 1])
text(5,50,'del')

subplot('position',[0.88 0.15 0.08 0.8])
ka=find(ismember(G.Type,'Amp'));
xa=[G.nlog10q(ka) G.nlog10q(ka)]'; xa=xa(:);
ya=[G.x1(ka) G.x2(ka)]'; ya=ya(:);
plot(xa,ya,'r','linewidth',2)
xlim([-.1 11])
set(gca,'ytick',xm(k),'yticklabel',[])
ylim([0 xe(end)])
set(gca,'ydir','reverse')
line([0;11]*ones(1,NC),[1;1]*xe,'color',0.5*[1 1 1],'linestyle','--')
text(0,-70,'Focal','fontsize',12,'horizontalAlign','center')
text(0,xe(end)+150,'-log10(q-value)','fontsize',12,'horizontalAlign','center')
line([0 11],[1 1],'linestyle','-','linewidth',1,'color',0*[1 1 1])
line([0 11 ],[1 1]*xe(end)-1,'linestyle','-','linewidth',1,'color',0*[1 1 1])
text(2.5,50,'gain')

% dx=0.5;
% for i=1:2:NC
%     s=-1
%     x1=[s-dx s+dx s+dx s-dx s-dx]';
%     y1=[xe(i) xe(i) xe(i+1) xe(i+1) xe(i)]'
%     patch(x1,y1,'k','edgecolor','none')   
% end

f=['~/Projects/Cancer/mouse/Gistic2/' SET '.SCNA.Fig2.' TODAY];
print(gcf,'-dpng','-r400',[f '.png'])
print(gcf,'-depsc','-r400',[f '.eps'])





