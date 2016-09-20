SIF=load_struct('~/amaro/deTiN_Figure_Code_and_PDFs/Figure4/SIF_for_cp_subclone.txt');

load('rb_colomap.mat')
load('grncolormap.mat')
figure()
hold on
starts=xhg19(1:22,0*(1:22),'hg19');
dx=.5;
for i=1:slength(SIF)
    s=i*2;
    acs=load_table(SIF.ACS_seg{i}); acs=rmfield(acs,{'header','headline'});
    deTiN=load_table(SIF.deTiN_segments{i}); deTiN=rmfield(deTiN,{'header','headline'});
    acs.t_n=acs.tau;
    acs.t_n(acs.tau>4)=4;
    acs.xStart=xhg19(acs.Chromosome,acs.Start_bp);
    acs.xEnd=xhg19(acs.Chromosome,acs.End_bp);
    
    if i==1
        for x=1:length(starts)
            if x<22
                l(x)=((starts(x+1)-starts(x))/2)+starts(x);
            else
                l(x)=starts(x)+((acs.xEnd(end)-starts(x))/2);
            end
        end

    end
    
    
    for j=1:slength(acs)
        y1=[ acs.xStart(j) acs.xStart(j) acs.xEnd(j) acs.xEnd(j) acs.xStart(j)];
        x1=[s-dx s+dx s+dx s-dx s-dx];
        c=rb_colormap(round((acs.t_n(j)/4)*255)+1,:);
        patch(x1,y1,c./255,'edgecolor','none')
    end
    z=s-1;
    for j=1:slength(deTiN)
        y1=[ deTiN.xs(j) deTiN.xs(j) deTiN.xe(j) deTiN.xe(j) deTiN.xe(j)];
        x1=[z-dx z+dx z+dx z-dx z-dx];
        k=find(acs.xStart==deTiN.xs(j));
        if deTiN.modal_TiN(j)>.5
            c=grncolor(end,:);
        else
        c=grncolor(round((deTiN.modal_TiN(j)/.5)*255)+1,:);
        end
        patch(x1,y1,c./255,'edgecolor','none')
        
    end
    
    plot([z-dx,z-dx],[0 l(end)+30000000],'k-')
    plot([s+dx,s+dx],[0 l(end)+30000000],'k-')
end

for x=1:length(starts)
            plot([0,slength(SIF)*2+1],[starts(x) starts(x)],'k--')
end
        set(gca,'ytick',l(1:2:22),'ytickLabel',num2str([1:2:22]'));
        ylim([-1 l(end)+30000000])
        box on
        xlim([.5 48.6])
        set(gca,'XTick',[])