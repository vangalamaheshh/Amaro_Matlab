function [percent_matching_segments, number_of_segments_seg1, number_of_segments_seg2...
    mean_seg_length]=compare_seg_files(segfile1,segfile2,outdir)

% load seg files in

seg_1=load_struct(segfile1);
seg_2=load_struct(segfile2);

%number of segments for each segment file

number_of_segments_seg1=slength(seg_1);
number_of_segments_seg2=slength(seg_2);


%if isequal(format1,'ReCapseg')  % no longer an issue due to column header
%fix
 %   [seg_1 seg_1_startx seg_1_endx]=matify_ReCapseg(seg_1); end

    [seg_1 seg_1_startx seg_1_endx]=matify_SNParray(seg_1); % makes columns vectors rather then cell arrays
%if isequal(format2,'ReCapseg')
%    [seg_2 seg_2_startx seg_2_endx]=matify_ReCapseg(seg_2); end

    [seg_2 seg_2_startx seg_2_endx]=matify_SNParray(seg_2); 
seg_1.ismatching=zeros(slength(seg_1),1);
seg_1.diff=nan(slength(seg_1),1);
for i=1:slength(seg_1)
    
    %%%%%%%%% check to see if any segment on the chromosome reciprocally overlaps 50%
    k=find(seg_2.Chromosome==seg_1.Chromosome(i));
     for j=1:length(k)
         sm=k(j);
         % make sure that the start is less than the end etc.
         if (seg_1.Start(i)<seg_2.End(sm) && seg_1.End(i)>seg_2.Start(sm))
     s=max([seg_1.Start(i);seg_2.Start(sm)]);
     e=min([seg_1.End(i);seg_2.End(sm)]);
     
     if ((e-s)/(seg_1.End(i)-seg_1.Start(i))) > .5 && ...
     ((e-s)/(seg_2.End(sm)-seg_2.Start(sm))) > .5
        seg_1.ismatching(i)=1;
        
        seg_1.diff(i)=seg_1.Segment_Mean(i)-seg_2.Segment_Mean(sm);
     
     end
         end
     end
         
    
end

figure()
% plotting both segment files along the genome
hold on
for i=1:length(seg_1.Segment_Mean)
line([seg_1_startx(i);seg_1_endx(i)],repmat(seg_1.Segment_Mean(i),2),'LineWidth',4,'Color',[255/255 99/255 71/255])
end

for i=1:length(seg_2.Segment_Mean)
line([seg_2_startx(i);seg_2_endx(i)],repmat(seg_2.Segment_Mean(i),2),'LineWidth',4,'Color',[71/255 99/255 255/255])
end

ylim([0 10])
strs=split(segfile1,'/');
s1name=strs{end};
strs=split(segfile2,'/');
s2name=strs{end};

text(mean(seg_1_startx),8,s1name,'Color',[255/255 99/255 71/255],'FontSize',12);
text(mean(seg_1_startx),7,s2name,'Color',[71/255 99/255 255/255],'FontSize',12);

saveas(gcf,strcat(outdir,s1name,'_vs_',s2name,'.fig'));
close(gcf)





seg_1.length=seg_1.End-seg_1.Start;
mean_matching_seg_length=mean(seg_1.length(seg_1.ismatching==1));
mean_seg_length=mean(seg_1.length);
%%%%%%%%%%%%%%%%%%%%%%%%%% replaced by unified break points method
%percent_matching_segments=sum(seg_1.ismatching)/slength(seg_1);
%
%agreement_between_matching_segments=sum(abs(seg_1.diff)<.1)/sum(seg_1.ismatching);
%old check for matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[uSeg1,uSeg2]=unify_break_points_recapseg(seg_1,seg_2);
% find difference between unified segment means
d=uSeg1.Segment_Mean-uSeg2.Segment_Mean;
% percent of unified segments that agree
percent_matching_segments=sum((abs(d)<.2))/slength(uSeg1);



end

function [seg_struct,xsegstart,xsegend]=matify_ReCapseg(seg_struct)
%% This function is no longer in use due to fixed columns in the RECASEG files

%     seg_struct.Chromosome=chrom2num(seg_struct.chrom);
 %    seg_struct.Start=str2double(seg_struct.locstart);
  %   seg_struct.End=str2double(seg_struct.locend);
   %  seg_struct.Segment_Mean=str2double(seg_struct.segmean);
    % for i=1:slength(seg_struct)
     %    seg_struct.Segment_Mean(i)=(2^seg_struct.Segment_Mean(i))*2;
     %end
     %xsegstart=xhg19(seg_struct.Chromosome,seg_struct.Start);
     %xsegend=xhg19(seg_struct.Chromosome,seg_struct.End);
end

function [seg_struct,xsegstart,xsegend]=matify_SNParray(seg_struct)
    seg_struct.Chromosome=chrom2num(seg_struct.Chromosome);
    seg_struct.Start=str2double(seg_struct.Start);
    seg_struct.End=str2double(seg_struct.End);
    seg_struct.Segment_Mean=str2double(seg_struct.Segment_Mean);
    for i=1:slength(seg_struct)
        seg_struct.Segment_Mean(i)=(2^seg_struct.Segment_Mean(i));
    end
    xsegstart=xhg19(seg_struct.Chromosome,seg_struct.Start);
    xsegend=xhg19(seg_struct.Chromosome,seg_struct.End);

end


function test
%segfile1='/Users/amaro/Downloads/segments.CESC-TCGA-FU-A23L-Tumor-SM-23G2X.tsv';
%segfile2='/Users/amaro/Downloads/LISTS_p_TCGA_129_134_148_149_S_GenomeWideSNP_6_H12_804844.nocnv_hg19.seg.txt';
%format1='ReCapseg';
%format2='SNParray';

files=load_struct('~/Projects/ReCapseg/STADOldCapsegReCapseg.txt');
for i=1:slength(files)
    segfile2=files.OldCapseg{i};
    format2='SNParray';
    segfile1=files.ReCapseg{i};
    format1='SNParray';
    [files.percent_matching_segments(i,1),files.nseg_recapseg(i,1),files.nseg_oldCapseg(i,1),files.prct_matching_that_agree(i,1),files.mean_RECAPSEG_length(i,1)]=compare_seg_files(segfile1,format1,segfile2,format2);
end

end
