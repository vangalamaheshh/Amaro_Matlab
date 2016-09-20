Table=load_struct('/Users/amaro/Downloads/ValidationHemangioma.txt');

N=slength(Table);
NV=3;



bottom_alignment=.77;
%% plotting purities
figure('position',[100 100 1000 1200])
hold on
subplot('position',[0.1 0.8 .89 0.15]) %left bottom width height
% P=bar([1:N],str2double([Table.purity_LOH	Table.purity_Del]),'stacked');
% set(P(2),'facecolor',[255/255 127/255 80/255],'EdgeColor',[1,1,1])
% set(P(1),'facecolor',[231/255 138/255 195/255],'EdgeColor',[1,1,1])
% set(gca,'xcolor',[1 1 1],'xtick',[]);
% ylabel('Purity'); box off
% Pbaseline=get(P,'BaseLine');
% set(Pbaseline{1},'Color',[.9 .9 .9],'LineWidth',1);
plot(str2double(Table.purity_Del),'x','MarkerFaceColor',[255/255 127/255 80/255],'MarkerEdgeColor',[255/255 127/255 80/255],'MarkerSize',14)
hold on
plot(str2double(Table.purity_LOH),'b.','MarkerSize',14)
set(gca,'xcolor',[1 1 1],'xtick',[]);
ylabel('Purity'); box off
hold on
subplot('position',[0.1 0.2 .895 0.54])
set(gca,'visible','off');
%axis([-1.5 N+2 -2 NV+4]);
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.1;


for v=1:NV
    for s=1:N
       
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])

        
    end
end


deletion_color=[43/255 131/255 186/255];
somatic_mut_color=[34/255 139/255 34/255];
Events={'VHL';'Del 3p';'LOH 3p';};



for v=1:NV
    for s=1:N
             x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];

       j=NV-(v-1);
        y_focal_del=[j-(dy_cn-.005) j-(dy_cn-.005) j j j-dy_cn];
        y_som=[j-(dy_cn-.005) j-(dy_cn-.005) j+(dy_cn-.005) j+(dy_cn-.005) j-dy_cn]; 
        
        y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn]; % arm lvl CN event
        
        if v==1;
        if isequal(Table.focalVHLdel{s},'1')
            patch(x1,y_focal_del,deletion_color,'edgecolor',deletion_color)
        end
        if isequal(Table.focalVHLdel{s},'2')
            patch(x1,y_focal_del,[.9 .9 .9],'edgecolor',deletion_color,'LineWidth',1)
        end
        if isequal(Table.VHL{s},'1')
            patch(x1,y_som,somatic_mut_color,'edgecolor',somatic_mut_color)
        end
        end
        
        if v==2;
        if isequal(Table.Del_chr3P{s},'1')
            patch(x1,y1,deletion_color,'edgecolor',deletion_color)
        end
        if isequal(Table.CN_neutral_LOH{s},'1')
            patch(x1,y1,[.9 .9 .9],'edgecolor',deletion_color,'LineWidth',1)
        end
        end
        
        if v==3;
            if isequal(Table.Homozygous{s},'2')
            patch(x1,y1,[178/255 171/255 210/255],'edgecolor',[178/255 171/255 210/255],'LineWidth',1)
            end
            if isequal(Table.Homozygous{s},'1')
            patch(x1,y_focal_del,[178/255 171/255 210/255],'edgecolor',[178/255 171/255 210/255],'LineWidth',1)
            end
        end
        
    end
end



