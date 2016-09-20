Recapseg=load_struct('/xchip/cga_home/amaro/CLL/ABSOLUTE/SCNAforcecall/Recapseg_segs.tsv');
SEGS=load_struct('/Users/amaro/Documents/SEGstoForceCall.txt');

for i=1:slength(Recapseg)
    Recapseg.SM{i,1}=Recapseg.sample_id{i}(end-4:end);
    Recapseg.Ind{i,1}=Recapseg.sample_id{i}(1:14);
end
for i=1:slength(SEGS)
    SEGS.sm{i,1}=SEGS.sample{i}(end-13:end-9);
    SEGS.individual_id{i,1}=SEGS.sample{i}(1:14);

end


Ind=unique(SEGS.individual_id);
for i=1:length(Ind)
    k=find(ismember(SEGS.individual_id,Ind{i}));
    kR=find(ismember(Recapseg.Ind,Ind{i}));
    clear s_f
    for j=1:length(kR)
        s=load_struct(Recapseg.recapseg_seg_file{kR(j)});
        s.Start=str2double(s.Start);
        s.End=str2double(s.End);
        for s_i=1:length(k)
        si=reorder_struct(s,ismember(s.Chromosome,SEGS.Chromosome{k(s_i)}));
        st=str2double(SEGS.Startbp{k(s_i)});
        ed=str2double(SEGS.Endbp{k(s_i)});
        if ~isempty(find(si.Start<=st&si.End>ed,1))
            k_i=find(si.Start<=st&si.End>=ed,1);
        else
            k_i=find(si.Start<=st,1,'last');
        end
    
        s_f.Sample{s_i,1}=Recapseg.sample_id{kR(j)};
        s_f.Chromsome{s_i,1}=SEGS.Chromosome{k(s_i)};
        s_f.Start{s_i,1}=SEGS.Startbp{k(s_i)};
        s_f.End{s_i,1}=SEGS.Endbp{k(s_i)};
        s_f.Num_Probes{s_i,1}=SEGS.n_probes{k(s_i)};
        s_f.Segment_Mean{s_i,1}=si.Segment_Mean{k_i};
       
        
        end
        save_struct(s_f,strcat('/xchip/cga_home/amaro/CLL/ABSOLUTE/SCNAforcecall/',s_f.Sample{1},'.forcecalled.seg'))
        clear s_f
    end

end
