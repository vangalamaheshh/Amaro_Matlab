%    output file columns:
%    1. line number in dranger results file
%    2. breakpoint number
%    3. first chromosome
%    4. predicted location on first chromosome
%    5. second chromosome
%    6. predicted location on second chromosome
%    7. is there an inversion?
%    8. how many reads support this breakpoints (in current runs stops after finding 3 supporting reads, so 3 means 3 or more)
%    9. the average score of the alignment of supporting reads
%   10. which sequence is first (0-first sequence 1-second sequence)
%   11. the length of micro homology
%   12. the length of foreign sequence added at the breakpoint
%   13. the foreign sequence added at the breakpoint

function estimate_breakpoint2(indlist,fileid,fastqname)

javaclasspath('/xchip/cga1/ydrier/findclusters/bwa/dist/sam-1.07.168.jar');
import net.sf.samtools.*;
import java.io.*;

outputf=fopen(['my_breakpoints_' fileid '.txt'],'w');

for idx=1:length(indlist)
    b=indlist(idx);
    if mod(idx,10)==1
        [idx b]
    end
    dirnum=num2str(b);
    fid=fopen([dirnum '/splitreads.helper']);
    t=fscanf(fid,'%d\t',9);
    seq1=fscanf(fid,'%s\n',1);
    seq2=fscanf(fid,'%s\n',1);
    fclose(fid);   
    inversion=t(4)==t(7);
    seq1=upper(seq1);
    seq2=upper(seq2);       
    chr=[t(2);t(5)];
    startongenome=[t(3);t(6)]-t(9);
    breakpointstart=[t(3);t(6)]-40;
    breakpointend=[t(3);t(6)]+50;
    if t(8)<=9
        breakpointstart=breakpointstart-40;
        breakpointend=breakpointend+340;
    end
    bprlen=breakpointend-breakpointstart+1;
    bpoffset=breakpointstart-startongenome;
    seq_start=max(bpoffset-101,1);
    seq_end=min(bpoffset+bprlen+101,([length(seq1); length(seq2)]));
    seq1=seq1(seq_start(1):seq_end(1));
    seq2=seq2(seq_start(2):seq_end(2));
    if inversion
        seq2=seqrcomplement(seq2);
    end

        [breakpoint,skip,score,seq2first,forseq,maxcounts]=...
  align_reads2(dirnum,fastqname,seq1,seq2,seq_start,seq_end,startongenome,t(4),t(7));

    if (score>0)
        overlap=0;
        j=breakpoint'-startongenome-seq_start+1;
        while ((j(1)-overlap>0)&&(j(1)-overlap<=length(seq1))&&(j(2)-overlap>0)&&(j(2)-overlap<=length(seq2))&&(seq1(j(1)-overlap)==seq2(j(2)-overlap)))
            overlap=overlap+1;
        end
        fprintf(outputf,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n',b,t(1),chr(1),breakpoint(1),chr(2),breakpoint(2),inversion,maxcounts,score,seq2first,overlap,skip,forseq);        
    end
end
fclose(outputf);


