SegmentFiles=dir('/xchip/cga_home/amaro/Mouse_SCLC/Segment_Files_Original_Set/');
MafFiles=dir('/xchip/cga_home/amaro/Mouse_SCLC/MAFS/');
SegmentFiles(1)=[];
SegmentFiles(1)=[];
MafFiles(1)=[];
MafFiles(1)=[];
Mycl1_Start=122671342;
Mycll_End=122681723;

for i=1:length(SegmentFiles)
 strs=regexp(SegmentFiles(i).name,'-','split');
SegmentFiles(i).uni_id=strs{6};
end

 for i=1:length(MafFiles)
for j=1:length(SegmentFiles)
if ~isempty(strfind(MafFiles(i).name,SegmentFiles(j).uni_id))
SegmentFiles(j).map=i;
end
end
end

counter=1;
counter2=1;
for i=1:length(SegmentFiles)
    seg_file=sprintf('/xchip/cga_home/amaro/Mouse_SCLC/Segment_Files_Original_Set/%s',SegmentFiles(i).name);
    maf_file=sprintf('/xchip/cga_home/amaro/Mouse_SCLC/MAFS/%s',MafFiles(SegmentFiles(i).map).name);
    Seg=load_struct(seg_file);
    Maf=load_struct(maf_file);
    Maf_no_4=trimStruct(Maf,~ismember(Maf.Chromosome,'4'));
    Maf_no_4=rmfields_if_exist(Maf,'N');
    MafChr4=trimStruct(Maf,ismember(Maf.Chromosome,'4'));
    MafChr4=rmfield_if_exist(MafChr4,'N');
    MafChr4.Start_position=cellfun(@str2num,MafChr4.Start_position);
    MafChr4.End_position=cellfun(@str2num,MafChr4.End_position);
    Seg=trimStruct(Seg,ismember(Seg.Chromosome,'chr4'));
    Seg.Start=cellfun(@str2num,Seg.Start);
    Seg.End=cellfun(@str2num,Seg.End);  
    Seg.Segment_Mean=cellfun(@str2num,Seg.Segment_Mean); 
    %Segment(1)=find(Seg.Start > Mycl1_Start,1,'last');
    Seg.Length=Seg.End-Seg.Start;
    Segment=find(Seg.Segment_Mean > .8 );
    
    mean(Seg.Segment_Mean)
%{ if Segment<length(Seg.Start) 
  %      if Seg.Segment_Mean(Segment+1)>.4
  %          Segment(3)=Segment(1)+1;
  %      end
  %  end
  %  if Seg.Segment_Mean(Segment-1)>.4
  %      Segment(2)=Segment(1)-1;
  %  end
    for k=1:length(Segment)
       Starts(k)=Seg.Start(Segment(k));
       Ends(k)=Seg.End(Segment(k));
       
    end
    Ending=max(Ends);
    Starting=min(Starts);
    mutations=find(MafChr4.Start_position>=Starting);
    check=find(MafChr4.End_position(mutations)<=Ending);
    mutations=mutations(check);
    length(MafChr4.Start_position)
    fields=fieldnames(Maf);
    length(mutations)
    for y=1:length(mutations)
        for x=1:length(fields)
            if iscell(MafChr4.(fields{x}))
        AllCHR4Mutations.(fields{x}){counter,1}=MafChr4.(fields{x}){mutations(y)};
            else
         AllCHR4Mutations.(fields{x})(counter,1)=MafChr4.(fields{x})(mutations(y));
            end

                
        end
        counter=counter+1;
    end
    
    i=1;
    for y=1:length(Maf.Chromosome)
        for x=1:length(fields)
            ALLOtherMutations.(fields{x}){counter2,1}=Maf_no_4.(fields{x}){y};
        end
        counter2=counter2+1;
    end
    
    
    
clear Segment
end

ALLOtherMutations.i_tumor_f=cellfun(@str2num,ALLOtherMutations.i_tumor_f);
AllCHR4Mutations.i_tumor_f=cellfun(@str2num,AllCHR4Mutations.i_tumor_f);
