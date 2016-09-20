
ABS_maf=load_struct('/Users/amaro/Documents/libfiltered_pon7_no_micro_sat_missing_removed.maf');
sig_genes=load_struct('/Users/amaro/Documents/sig_genes_Louis.txt');

 missense=[51,160,44]./255; %green
 nonsense= [227,26,28]./255; %red
 splice= [255,127,0]./255; %orange
 indel= [106,61,154]./255; %purple
 sig_genes.q=str2double(sig_genes.q);
 sig_genes=reorder_struct(sig_genes,sig_genes.q<.1);
samples.pair_id=unique(ABS_maf.Tumor_Sample_Barcode);

 ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Hugo_Symbol,sig_genes.gene));
 ABS_maf=reorder_struct(ABS_maf,ismember(ABS_maf.Variant_Classification,'Missense_Mutation')|ismember(ABS_maf.Variant_Classification,'Nonsense_Mutation')|...
     ismember(ABS_maf.Variant_Classification,'Splice_Site')|ismember(ABS_maf.Variant_Classification,'Frame_Shift_Ins')...
     |ismember(ABS_maf.Variant_Classification,'Frame_Shift_Del')|ismember(ABS_maf.Variant_Classification,'In_Frame_Del')...
     |ismember(ABS_maf.Variant_Classification,'In_Frame_Ins'));
 
 N=slength(samples);
 NV=slength(sig_genes);
 subplot('position',[0.1 0.05 0.85 0.7]) %left bottom width hieght
axis([-5 N+3 -2 NV+3]);
set(gca,'visible','off');
%setting dims of boxs for amps and deletions
dx_cn=0.4; dy_cn=0.4;
dx_mut=0.4; dy_mut=.4;

 disp(strcat('Plotting'));

for i=1:slength(sig_genes)
    
   
    for j=1:slength(samples)
        
        k=ismember(ABS_maf.Tumor_Sample_Barcode,samples.pair_id{j});
        if ismember(sig_genes.gene{i},ABS_maf.Hugo_Symbol(k))
                    samples.(sig_genes.gene{i})(j,1)=1;
        else
        samples.(sig_genes.gene{i})(j,1)=0;
        end
    end
end


f_samples=fieldnames(samples);
samples=sort_struct(samples,{f_samples{2:end}},ones(length(f_samples)-1,1)*-1);


for v=1:NV
    for s=1:N
       
     x1=[s-dx_cn s+dx_cn s+dx_cn s-dx_cn s-dx_cn];
    j=NV-(v-1);
    y1=[j-dy_cn j-dy_cn j+dy_cn j+dy_cn j-dy_cn];
   
    patch(x1,y1,[.9 .9 .9],'edgecolor',[.9 .9 .9])

        
    end
end

for v=1:NV
    for s=1:slength(samples)
    
        
    x1=[s-dx_mut s+dx_mut s+dx_mut s-dx_mut s-dx_mut];
    j=NV-(v-1);
    y1=[j-dy_mut j-dy_mut j+dy_mut j+dy_mut j-dy_mut];
    
    k=ismember(ABS_maf.Tumor_Sample_Barcode,samples.pair_id{s});  
    if ismember(sig_genes.gene{v},ABS_maf.Hugo_Symbol(k))
         l=find(ismember(ABS_maf.Hugo_Symbol,sig_genes.gene{v})&ismember(ABS_maf.Tumor_Sample_Barcode,samples.pair_id{s}));
          
    if isequal(ABS_maf.Variant_Type{l},'DEL')||isequal(ABS_maf.Variant_Type{l},'INS')
          patch(x1,y1,indel,'edgecolor',indel)
    elseif isequal(ABS_maf.Variant_Classification{l},'Missense_Mutation')
                  patch(x1,y1,missense,'edgecolor',missense)
    elseif isequal(ABS_maf.Variant_Classification{l},'Nonsense_Mutation')
                  patch(x1,y1,nonsense,'edgecolor',nonsense)
    elseif isequal(ABS_maf.Variant_Classification{l},'Splice_Site')
                  patch(x1,y1,splice,'edgecolor',splice)
    end
                    ss(v,s)=1;
                    
    end
    end
end
    
%labels
for v=1:NV
text (-.5,NV-(v-.4),sig_genes.gene{v},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16)
end


