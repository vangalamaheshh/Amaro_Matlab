function landscape_for_germ_maf(germline_maf,selected_genes,output_file)
set(0,'DefaultFigureVisible','off')
disp('Loading Maf');
G=load_struct(germline_maf);
s=load_struct_noheader(selected_genes);

si.gene=s.col1;

G=reorder_struct(G,ismember(G.gene,si.gene));
[tr si.gene]=count(G.gene,-1);
%plotting

N=length(unique(G.sample));
NV=slength(si);

samples.names=unique(G.sample);

disp('Sorting');
% sorting code
for v=1:NV
    for s=1:N
         k=ismember(G.sample,samples.names{s});
         if ismember(si.gene{v},G.gene(k));
             samples.(si.gene{v})(s,1)=1;
         else
             samples.(si.gene{v})(s,1)=0;
         end
        
    end
end
fields=fieldnames(samples);
samples=sort_struct(samples,{fields{2:end}},(ones(length({fields{2:end}}))*-1));


% current set up to match oncoprint
somatic_mut_color=[34/255 139/255 34/255];
germ_mut_color=[138/255 43/255 226/255];
deletion_color=[0 0 1];



hold on


N=length(unique(G.sample));
NV=slength(si);

axis([-1.5 N+2 -2 NV+4]);
set(gca,'visible','off');
%setting dims of boxs 
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.33; dy_mut=.35;
disp('Plotting');
% place base line patches
for v=1:NV
    for s=1:N
        
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])

        
    end
end


% place germline / somatic / LOH symbols
for v=1:NV
    for s=1:N
        
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];%
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
    
    y_germ=[j-(dy_cn-.005) j-(dy_cn-.005) j j j-dy_cn];% germline y
    y_som=[j j j+(dy_cn-.005) j+(dy_cn-.005) j]; %somatic y
    
    k=ismember(G.sample,samples.names{s});
    if ismember(si.gene{v},G.gene(k));
        l=find(ismember(G.gene,si.gene{v})&ismember(G.sample,samples.names{s}));
        for i=1:length(l)
            if isequal(G.loh{l(i)},'True')||isequal(G.loh{l(i)},'TRUE')
            patch(x1,y1,[.9 .9 .9],'edgecolor',deletion_color,'LineWidth',2)

           
               
            end
            
            if isequal(G.is_somatic{l(i)},'TRUE')|| isequal(G.is_somatic{l(i)},'True')
                patch(x1,y_som,somatic_mut_color,'edgecolor',somatic_mut_color,'EdgeAlpha',0)
            else
                patch(x1,y_germ,germ_mut_color,'edgecolor',germ_mut_color,'EdgeAlpha',0)
            end
        
        end
            
        
    end
        

    end
end


for v=1:NV
text (0,NV-(v-.65),si.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end

print_D(output_file,{{'png','-r300'}});
close all 

end

function test
germline_maf='/Users/amaro/Downloads/LandScapes/example.maf';
selected_genes='/Users/amaro/Downloads/LandScapes/geneset.txt';

landscape_for_germ_maf(germline_maf,selected_genes,'output');


end