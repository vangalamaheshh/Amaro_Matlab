external_ids=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/full_list_of_external_samples.tsv');
%SEG=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/AggregatedCCFCLL8.called.seg'); 
SEG=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/AggregatedCCFCLL8.called.3.18.seg');
pair_map=load_struct('/Volumes/xchip_cga_home/amaro/CLL/Stilgenbauer.DFCI.ICGC.630.pairs.tsv');
external_ids=reorder_struct(external_ids,~ismember(external_ids.external_id_capture,''));
overlap = @(x1,x2,y1,y2) max([0, min([x2;y2]) - max([x1;y1])]);
seg_pairs=unique(SEG.sample);
for i=1:length(seg_pairs)
    k=ismember(SEG.sample,seg_pairs{i});
    l=ismember(pair_map.pair_id,seg_pairs{i});
    SEG.sample_id(k,1)=pair_map.case_sample(l);
    SEG.individual(k,1)=pair_map.patient_id(l);
end



for i=1:slength(external_ids)
external_ids.TP{i,1}=external_ids.external_id_capture{i}(end-1:end);
external_ids.sm_id{i,1}=external_ids.sample_id{i}(end-7:end);
end

external_ids=reorder_struct(external_ids,ismember(external_ids.TP,{'01','02'}));
SEG=reorder_struct(SEG,ismember(SEG.sample_id,external_ids.sample_id));
cheatsheet=load_struct('/Users/amaro/Downloads/c.txt');

for i=1:slength(SEG)
    SEG.individual{i}=SEG.sample{i}(5:13);
end



SEG=reorder_struct(SEG,ismember(SEG.individual,cheatsheet.sample));
events=unique(SEG.ID);
events=events(~ismember(events,'NaN'));
samples=unique(SEG.sample_id);

for i=1:length(samples)
    k=ismember(SEG.sample_id,samples{i});
    l=ismember(external_ids.sample_id,samples{i});
    SEG.TP(k,1)=external_ids.TP(l);
end

individuals=unique(SEG.individual);
SEG.CCF=str2double(SEG.CCF);
SEG.CCF_High=str2double(SEG.CCF_High);
SEG.CCF_Low=str2double(SEG.CCF_Low);
SEG.Start_bp=str2double(SEG.Start_bp);
SEG.End_bp=str2double(SEG.End_bp);
SEG.xStart=xhg19(SEG.Chromosome,SEG.Start_bp);
SEG.xEnd=xhg19(SEG.Chromosome,SEG.End_bp);



for j=1:length(individuals)
    i_SEG=reorder_struct(SEG,ismember(SEG.individual,individuals{j}));
    for i=1:length(events)
        
            
            
            %%%%%%%%%%% check TP1
            if ismember(events{i},i_SEG.ID(ismember(i_SEG.TP,'01')))
                SCNA_table.(strcat(events{i},'TP1'))(j,1)=max(i_SEG.CCF(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'01')));
                SCNA_table.(strcat(events{i},'TP1_High'))(j,1)=max(i_SEG.CCF_High(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'01')));
                SCNA_table.(strcat(events{i},'TP1_Low'))(j,1)=max(i_SEG.CCF_Low(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'01')));
            elseif ismember(events{i},i_SEG.ID(ismember(i_SEG.TP,'02')))
                
                [l v]=max(i_SEG.CCF(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'02')));
                starts=i_SEG.xStart(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'02'));
                ends=i_SEG.xEnd(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'02'));
                max_ov=0;
                for z=1:slength(i_SEG)
                    
                    if overlap(starts(v),ends(v),i_SEG.xStart(z),i_SEG.xEnd(z))>max_ov && isequal('01',i_SEG.TP{z})
                        max_ov=overlap(starts(v),ends(v),i_SEG.xStart(z),i_SEG.xEnd(z));
                        max_ov_k=z;
                    end
                end
                
                SCNA_table.(strcat(events{i},'TP1'))(j,1)=i_SEG.CCF(max_ov_k);
                SCNA_table.(strcat(events{i},'TP1_Low'))(j,1)=i_SEG.CCF_Low(max_ov_k);
                SCNA_table.(strcat(events{i},'TP1_High'))(j,1)=i_SEG.CCF_High(max_ov_k);
            end
            
            %%%%%%%%% Check TP2
            if ismember(events{i},i_SEG.ID(ismember(i_SEG.TP,'02')))
                SCNA_table.(strcat(events{i},'TP2'))(j,1)=max(i_SEG.CCF(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'02')));
                SCNA_table.(strcat(events{i},'TP2_Low'))(j,1)=max(i_SEG.CCF_Low(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'02')));
                SCNA_table.(strcat(events{i},'TP2_High'))(j,1)=max(i_SEG.CCF_High(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'02')));
                
                
            elseif ismember(events{i},i_SEG.ID(ismember(i_SEG.TP,'01')))
                [l v]=max(i_SEG.CCF(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'01')));
                starts=i_SEG.xStart(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'01'));
                ends=i_SEG.xEnd(ismember(i_SEG.ID,events{i})&ismember(i_SEG.TP,'01'));
                max_ov=0;
                for z=1:slength(i_SEG)
                    
                    if overlap(starts(v),ends(v),i_SEG.xStart(z),i_SEG.xEnd(z))>max_ov && isequal('02',i_SEG.TP{z})
                        max_ov=overlap(starts(v),ends(v),i_SEG.xStart(z),i_SEG.xEnd(z));
                        max_ov_k=z;
                    end
                end
                SCNA_table.(strcat(events{i},'TP2'))(j,1)=i_SEG.CCF(max_ov_k);
                SCNA_table.(strcat(events{i},'TP2_Low'))(j,1)=i_SEG.CCF_Low(max_ov_k);
                SCNA_table.(strcat(events{i},'TP2_High'))(j,1)=i_SEG.CCF_High(max_ov_k);
            end
            
            %%%%%%%% Check if not in I_Seg at all
            if ~ismember(events{i},i_SEG.ID(ismember(i_SEG.TP,'02'))) && ~ismember(events{i},i_SEG.ID(ismember(i_SEG.TP,'01')))
                SCNA_table.(strcat(events{i},'TP1'))(j,1)=NaN;
                SCNA_table.(strcat(events{i},'TP2'))(j,1)=NaN;
                SCNA_table.(strcat(events{i},'TP1_Low'))(j,1)=NaN;
                SCNA_table.(strcat(events{i},'TP2_Low'))(j,1)=NaN;
                SCNA_table.(strcat(events{i},'TP1_High'))(j,1)=NaN;
                SCNA_table.(strcat(events{i},'TP2_High'))(j,1)=NaN;
            end
            
            
        end
   

end


for i=1:length(events)
    hold on
    num_val=sum(~isnan(SCNA_table.(strcat(events{i},'TP1'))));
    tp1x=.8:.4/num_val:1.2;
    tp2x=1.8:.4/num_val:2.2;
    count=1;
    SCNA_table=sort_struct(SCNA_table,(strcat(events{i},'TP1')),-1);
    SCNA_table=sort_struct(SCNA_table,(strcat(events{i},'TP2')),-1);
    for j=1:length(individuals)
        if ~isnan(SCNA_table.(strcat(events{i},'TP1'))(j))
            if SCNA_table.(strcat(events{i},'TP1_High'))(j)<SCNA_table.(strcat(events{i},'TP2_Low'))(j) %% Red
                errorbar([tp1x(count),tp2x(count)],[SCNA_table.(strcat(events{i},'TP1'))(j),SCNA_table.(strcat(events{i},'TP2'))(j)],...
                    [(SCNA_table.(strcat(events{i},'TP1'))(j)-SCNA_table.(strcat(events{i},'TP1_Low'))(j)),...
                    (SCNA_table.(strcat(events{i},'TP2'))(j)-SCNA_table.(strcat(events{i},'TP2_Low'))(j))],...
                    [(SCNA_table.(strcat(events{i},'TP1_High'))(j)-SCNA_table.(strcat(events{i},'TP1'))(j)),...
                    (SCNA_table.(strcat(events{i},'TP2_High'))(j)-SCNA_table.(strcat(events{i},'TP2'))(j))],'r')
                
            elseif SCNA_table.(strcat(events{i},'TP1_Low'))(j)>SCNA_table.(strcat(events{i},'TP2_High'))(j) %% Blue
                errorbar([tp1x(count),tp2x(count)],[SCNA_table.(strcat(events{i},'TP1'))(j),SCNA_table.(strcat(events{i},'TP2'))(j)],...
                    [(SCNA_table.(strcat(events{i},'TP1'))(j)-SCNA_table.(strcat(events{i},'TP1_Low'))(j)),...
                    (SCNA_table.(strcat(events{i},'TP2'))(j)-SCNA_table.(strcat(events{i},'TP2_Low'))(j))],...
                    [(SCNA_table.(strcat(events{i},'TP1_High'))(j)-SCNA_table.(strcat(events{i},'TP1'))(j)),...
                    (SCNA_table.(strcat(events{i},'TP2_High'))(j)-SCNA_table.(strcat(events{i},'TP2'))(j))],'b')
            else %Grey
                errorbar([tp1x(count),tp2x(count)],[SCNA_table.(strcat(events{i},'TP1'))(j),SCNA_table.(strcat(events{i},'TP2'))(j)],...
                    [(SCNA_table.(strcat(events{i},'TP1'))(j)-SCNA_table.(strcat(events{i},'TP1_Low'))(j)),...
                    (SCNA_table.(strcat(events{i},'TP2'))(j)-SCNA_table.(strcat(events{i},'TP2_Low'))(j))],...
                    [(SCNA_table.(strcat(events{i},'TP1_High'))(j)-SCNA_table.(strcat(events{i},'TP1'))(j)),...
                    (SCNA_table.(strcat(events{i},'TP2_High'))(j)-SCNA_table.(strcat(events{i},'TP2'))(j))],'Color',[.5 .5 .5])
                
            end
            count=count+1;
        end
    end
    title(events{i},'FontSize',28)
    xlim([0.5 2.5])
    ylim([0 1]);
    ylabel('Cancer Cell Fraction','FontSize',20)
    set(gca,'YTick',[0;.5;1],'XTick',[1;2],'XTickLabel',{'TP1';'TP2'})
    print(gcf,'-depsc',strcat('/Users/amaro/Documents/GeneProgressionPlotsCLLStilgenbauer/CNVs/',events{i},'.eps'));
    close all
end


