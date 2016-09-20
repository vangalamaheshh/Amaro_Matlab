amps=load_struct('GCT_Armlvl_data/amp_list_march.txt');
dels=load_struct('GCT_Armlvl_data/del_list_march.txt');
cnloh=load_struct('GCT_Armlvl_data/cnloh_list_march.txt');
arm_locs=load_struct_noheader('GCT_Armlvl_data/cytoBandArms.txt');


starts=xhg19(chromosome2num_legacy(arm_locs.col1),str2double(arm_locs.col2));
starts=sort(starts);
ends=xhg19(chromosome2num_legacy(arm_locs.col1),str2double(arm_locs.col3));
ends=sort(ends);
starts(end+1)=ends(end);
ff=fieldnames(amps);

ClinicalData=load_struct('GCT_Armlvl_data/ClinicalPathData.txt');
for i=1:slength(ClinicalData)
ClinicalData.sample{i,1}=regexp(ClinicalData.FH_names{i},'DFCI_[0-9]+','match');
ClinicalData.sample{i,1}=char(ClinicalData.sample{i});
end
for i=2:length(ff)
arms.sample{i-1,1}=regexp(ff{i},'DFCI_[0-9]+','match');
arms.sample{i-1,1}=char(arms.sample{i-1,1});
end
amps=rename_fields(amps,{ff{2:end}},arms.sample);
dels=rename_fields(dels,{ff{2:end}},arms.sample);
cnloh=rename_fields(cnloh,{ff{2:end}},arms.sample);

ClinicalData=reorder_struct(ClinicalData,ismember(ClinicalData.sample,arms.sample));

figure()
hold on
%xlim([0.65 41])
set(gca,'ytick',starts(2:2:37),'ytickLabel',arm_locs.col4([2:2:slength(arm_locs)-2]'));
set(gca,'xtick',[])
%xlabel('Samples','FontSize',20)
ylabel('Chr','FontSize',20)
%ylim([0 starts(end-2)])
for j=1:length(starts)-1
    i=1;
    dx=.3;
    x1=[i-dx i+dx i+dx i-dx i-dx];
    y1=[starts(j) starts(j) starts(j+1) starts(j+1) starts(j)];
    if mod(j,2)==0
        patch(x1,y1,[0 0 0],'edgecolor','none')
    else
        patch(x1,y1,[0.7 0.7 0.7],'edgecolor','none')
    end
    
    for i=1:slength(ClinicalData)
        
        dx=.5;
        
        
        x1=[1+i-dx 1+i+dx 1+i+dx 1+i-dx 1+i-dx];
        y1=[starts(j) starts(j) starts(j+1) starts(j+1) starts(j)];
        
        if str2double(amps.(ClinicalData.sample{i})(j))==1
          %  patch(x1,y1,[1 0 0],'edgecolor','none')
        end
        
        if str2double(dels.(ClinicalData.sample{i})(j))==1
         %  patch(x1,y1,[0 0 1],'edgecolor','none')
        end
        
        if str2double(cnloh.(ClinicalData.sample{i})(j))==1
          patch(x1,y1,[50/255,205/255,50/255],'edgecolor','none')
        end
        
        
        
        
    end
end
extra_track_start=starts(end);
for i=1:slength(ClinicalData)
    s=i+1+.5;
    x1=[s-dx s+dx s+dx s-dx s-dx];
    y1=[ extra_track_start extra_track_start extra_track_start+100000000 extra_track_start+100000000 extra_track_start];
    patch(x1,y1,[.5 .5 .5],'edgecolor','none')
    if isequal(ClinicalData.Vitalstatus{i},'1')
        patch(x1,y1,[1 0 0],'edgecolor','none')
    else
        patch(x1,y1,[0 0 0],'edgecolor','none')
    end
end

xticklabel_rotate_simplified([1:slength(ClinicalData)]+.5,90,ClinicalData.sample)

