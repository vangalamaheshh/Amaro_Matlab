%compare break points of most related and most unrelated pairs

overlap = @(x1,x2,y1,y2) max([0, min([x2;y2]) - max([x1;y1])]);
% 
% 
% SEG=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/ACS_frozen.seg');
% SEG.xs=xhg19(chromosome2num_legacy(SEG.Chromosome),str2double(SEG.Startbp));
% SEG.xe=xhg19(chromosome2num_legacy(SEG.Chromosome),str2double(SEG.Endbp));
% SEG.tau=str2double(SEG.tau);
% SEG.n_hets=str2double(SEG.n_hets);
SampleTable=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/SamplesBarretsFrozenTrios.tsv');
SampleTable=reorder_struct(SampleTable,ismember(SampleTable.sample_type,'Normal'));

ABS_seg_file=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/Aggregate_Seg_file.txt');
ABS_seg_file=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/AggregatedFFPE_abs.seg');

targets=load_table('/xchip/cga_home/amaro/Barretts/ABSOLUTE_Run_Paper_Combined_Sets/BreakpointAnalysis/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets_germline_copy_number_variants_X_Y_removed.bed');
ABS_seg_file.corrected_total_cn=str2double(ABS_seg_file.corrected_total_cn);
targets.xs=xhg19(targets.chr,targets.start);
targets.xe=xhg19(targets.chr,targets.end);
targets.xmid=targets.xs+((targets.xe-targets.xs)/2);
ABS_seg_file.xs=xhg19(ABS_seg_file.Chromosome,str2double(ABS_seg_file.Startbp));
ABS_seg_file.xe=xhg19(ABS_seg_file.Chromosome,str2double(ABS_seg_file.Endbp));


%samples=load_struct('/Users/amaro/Documents/BarrettsRevisionFigures/TSPs_Sample_Map.txt');
samples=load_struct('Users/amaro/Downloads/BarrettsFFPECaseMappingFile.txt')
samples.relatedness=str2double(samples.relatedness);
unrelated_ind=unique(samples.individual_id(samples.relatedness==0));
related_ind=unique(samples.individual_id(samples.relatedness>.3));
individuals=unique(samples.individual_id);
samples.shared_segments=zeros(slength(samples),1);
samples.shared_breakpoints=zeros(slength(samples),1);


for j=1:length(individuals)
clear b_seg e_seg
b_seg=reorder_struct(ABS_seg_file,ismember(ABS_seg_file.sample,samples.pair_id(ismember(samples.individual_id,individuals{j})&~ismember(samples.tissue,'ESO'))));
e_seg=reorder_struct(ABS_seg_file,ismember(ABS_seg_file.sample,samples.pair_id(ismember(samples.individual_id,individuals{j})&ismember(samples.tissue,'ESO'))));
n_seg=load_struct(SampleTable.recapseg_seg_file{ismember(SampleTable.individual_id,individuals{j})});
n_seg.xb=xhg19(n_seg.Chromosome,str2double(n_seg.Start));
[breakpoints_e, breakpoints_b]=segs_to_breakpoints(e_seg,b_seg,targets,.5,n_seg);


shared=0;

for i=1:slength(breakpoints_e)
    if ~isempty(find(breakpoints_b.target_range_start<breakpoints_e.pos(i) & breakpoints_e.pos(i)<breakpoints_b.target_range_end))
        
        shared=shared+1;
        shared_e_br.xp(shared,1)=breakpoints_e.pos(i);
        shared_e_br.Chr(shared,1)=e_seg.Chromosome(ismember(e_seg.xs,breakpoints_e.pos(i))|ismember(e_seg.xe,breakpoints_e.pos(i)));
        
        if ~isempty(e_seg.Startbp(ismember(e_seg.xs,breakpoints_e.pos(i))))
            shared_e_br.pos(shared,1)=e_seg.Startbp(ismember(e_seg.xs,breakpoints_e.pos(i)));
        else
            shared_e_br.pos(shared,1)=e_seg.Endbp(ismember(e_seg.xe,breakpoints_e.pos(i)));
        end
    end
end
    
samples.shared_breakpoints(ismember(samples.individual_id,individuals{j})&ismember(samples.tissue,'ESO'),1)=shared;
shared=0;

for i=1:slength(breakpoints_b)
    if ~isempty(find(breakpoints_e.target_range_start<breakpoints_b.pos(i) & breakpoints_b.pos(i)<breakpoints_e.target_range_end))
        shared=shared+1;
        shared_b_br.xp(shared,1)=breakpoints_b.pos(i);
        shared_b_br.Chr(shared,1)=b_seg.Chromosome(ismember(b_seg.xs,breakpoints_b.pos(i))|ismember(b_seg.xe,breakpoints_b.pos(i)));
        
        if ~isempty(b_seg.Startbp(ismember(b_seg.xs,breakpoints_b.pos(i))))
            shared_b_br.pos(shared,1)=b_seg.Startbp(ismember(b_seg.xs,breakpoints_b.pos(i)));
        else
            shared_b_br.pos(shared,1)=b_seg.Endbp(ismember(b_seg.xe,breakpoints_b.pos(i)));
        end
        
    end
end

samples.shared_breakpoints(ismember(samples.individual_id,individuals{j})&~ismember(samples.tissue,'ESO'),1)=shared;

samples.total_breakpoints(ismember(samples.individual_id,individuals{j})&ismember(samples.tissue,'ESO'),1)=slength(breakpoints_e);
samples.total_breakpoints(ismember(samples.individual_id,individuals{j})&~ismember(samples.tissue,'ESO'),1)=slength(breakpoints_b);



breakpoint_matrix.chromosome=breakpoints_e.Chromosome;
breakpoint_matrix.position=breakpoints_e.Breakpoint_Position;
breakpoint_matrix.target_range_start=breakpoints_e.target_range_start;
breakpoint_matrix.target_range_end=breakpoints_e.target_range_end;
breakpoint_matrix.EAC=ones(slength(breakpoint_matrix),1);
breakpoint_matrix.BES=zeros(slength(breakpoint_matrix),1);
breakpoint_matrix.EAC_dCN=breakpoints_e.dCN;
breakpoint_matrix.BE_dCN=zeros(slength(breakpoint_matrix),1);


num_bps=slength(breakpoint_matrix);

for i=1:slength(breakpoints_b)
    if isempty(find(breakpoints_e.target_range_start<breakpoints_b.pos(i) & breakpoints_b.pos(i)<breakpoints_e.target_range_end))
        num_bps=num_bps+1;
        breakpoint_matrix.chromosome{num_bps}=breakpoints_b.Chromosome{i};
        breakpoint_matrix.position{num_bps}=breakpoints_b.Breakpoint_Position{i};
        breakpoint_matrix.target_range_start(num_bps)=breakpoints_b.target_range_start(i);
        breakpoint_matrix.target_range_end(num_bps)=breakpoints_b.target_range_end(i);
        breakpoint_matrix.BES(num_bps)=1;
        breakpoint_matrix.EAC(num_bps)=0;
        breakpoint_matrix.BE_dCN(num_bps)=breakpoints_b.dCN(i);
        breakpoint_matrix.EAC_dCN(num_bps)=0;
    else
        breakpoint_matrix.BES(breakpoints_e.target_range_start<breakpoints_b.pos(i) & breakpoints_b.pos(i)<breakpoints_e.target_range_end)=1;
        breakpoint_matrix.BE_dCN(breakpoints_e.target_range_start<breakpoints_b.pos(i) & breakpoints_b.pos(i)<breakpoints_e.target_range_end)=breakpoints_b.dCN(i);
    end    
end

breakpoint_matrix=sort_struct(breakpoint_matrix,'target_range_start');

save_struct(breakpoint_matrix,strcat('/Users/amaro/Documents/BarrettsRevisionFigures/',samples.individual_id{ismember(samples.individual_id,individuals{j})&~ismember(samples.tissue,'ESO')},'_breakpoints_matrix.txt'));
breakpoints_b=rmfield(breakpoints_b,'pos');
breakpoints_e=rmfield(breakpoints_e,'pos');
save_struct(breakpoints_b,strcat('/Users/amaro/Documents/BarrettsRevisionFigures/',samples.pair_id{ismember(samples.individual_id,individuals{j})&~ismember(samples.tissue,'ESO')},'_breakpoints.txt'));
save_struct(breakpoints_e,strcat('/Users/amaro/Documents/BarrettsRevisionFigures/',samples.pair_id{ismember(samples.individual_id,individuals{j})&ismember(samples.tissue,'ESO')},'_breakpoints.txt'));

clear breakpoint_matrix



end









% 
% 
% 
% 
% for i=1:length(individuals)
% 
% b_seg=reorder_struct(SEG,ismember(SEG.Sample,samples.pair_id(ismember(samples.individual_id,individuals{i})&~ismember(samples.tissue,'ESO'))));
% e_seg=reorder_struct(SEG,ismember(SEG.Sample,samples.pair_id(ismember(samples.individual_id,individuals{i})&ismember(samples.tissue,'ESO'))));
% for j=1:slength(e_seg)
%     if e_seg.tau(j)<=1.5 || e_seg.tau(j)>=2.5
%     [distance_s, index_s]=min(abs(b_seg.xs-e_seg.xs(j)));
%     [distance_e, index_e]=min(abs(b_seg.xe-e_seg.xe(j)));
%     
%     if distance_s<distance_e
%         index=index_s;
%     else
%         index=index_e;
%     end
%     
%     s_diff=abs(b_seg.xs(index)-e_seg.xs(j));
%     e_diff=abs(b_seg.xe(index)-e_seg.xe(j));
%     segment_overlap=overlap(b_seg.xs(index),b_seg.xe(index),e_seg.xs(j),e_seg.xe(j));
%     overlap_b=segment_overlap/(b_seg.xe(index)-b_seg.xs(index));
%     overlap_e=segment_overlap/(e_seg.xe(j)-e_seg.xs(j));
%     if overlap_b>=.51 && overlap_e>=.51 && ((b_seg.tau(index)<=1.5 && e_seg.tau(j)<=1.5)||(b_seg.tau(index)>=2.5&&e_seg.tau(j)>=2.5)) && e_seg.n_hets(j)>3
%         samples.shared_segments(ismember(samples.individual_id,individuals{i}))=samples.shared_segments(ismember(samples.individual_id,individuals{i}))+1;
%     end
%     end
% 
% 
% 
% end
% end

