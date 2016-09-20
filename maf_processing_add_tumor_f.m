case_control=load_struct('/xchip/cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/case_control_table.tsv');
external_id=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/external_id_capture');

for i=1:slength(external_id)
 k=ismember(IGVmaf.Matched_Norm_Sample_Barcode,external_id.external_id_capture{i});    
if sum(k)>0
p_l=find(ismember(case_control.case_sample,external_id.sample_id{i}));
if isempty(p_l)
p_l=find(ismember(case_control.control_sample,external_id.sample_id{i}));
end
if ~isempty(p_l)
IGVmaf.Matched_Norm_Sample_Barcode(k)=case_control.control_sample(p_l);
end
end
end



SNPs=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/StilgenbauerSNPs.clean.maf');
SNPs=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/MSNP_ICGC.aggregate.maf');



SNPs.key=strcat(SNPs.Chromosome,SNPs.Start_position,SNPs.Tumor_Sample_Barcode);

[i m]=ismember(SNPs.key,IGVmaf.key);                                                                                                                                          
 IGVmaf.i_tumor_f(m(m>0),1)=SNPs.i_tumor_f(i);  
 IGVmaf.t_alt_count(m(m>0),1)=SNPs.t_alt_count(i);
 IGVmaf.t_ref_count(m(m>0),1)=SNPs.t_ref_count(i); 
 
 M2=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/StilgenbauerM2.clean.maf'); 
 M2=load_struct('~/Projects/CLL/StilgenbauerMafFreeze2.0/M2_ICGC.aggregate.maf');
 
 
 for i=1:slength(external_id)
k=ismember(M2.Matched_Norm_Sample_Barcode,external_id.external_id_capture{i});

if sum(k)>0     
p_l=find(ismember(case_control.control_sample,external_id.sample_id{i}),1);
if isempty(p_l)
p_l=find(ismember(case_control.case_sample,external_id.sample_id{i}),1);   
end
if ~isempty(p_l)
M2.Tumor_Sample_Barcode(k)=case_control.case_sample(p_l);              
end
end
end

 
 
 M2.key=strcat(M2.Chromosome,M2.Start_position,M2.Tumor_Sample_Barcode);        
 
 
 for i=1:slength(M2)
str=split(M2.allelic_depth{i},',');
M2.t_ref_count{i,1}=str{1};
M2.t_alt_count{i,1}=str{2};
 end
[i m]=ismember(M2.key,IGVmaf.key); 

  IGVmaf.t_alt_count(m(m>0),1)=M2.t_alt_count(i); 
 IGVmaf.t_ref_count(m(m>0),1)=M2.t_ref_count(i); 
 
 
 
 
 
 
 
 
 
 IGVmaf=reorder_struct(IGVmaf,~ismember(IGVmaf.Tumor_Sample_Barcode,'CLL-GCLL-0112-Tumor-SM-41JXE'));
 