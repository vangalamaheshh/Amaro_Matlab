function D = getdict(bamfile)
% getdict(bamfile)
%
% returns a cell array of strings holding the sequence names
% of the specified BAM file's dictionary.
%
% position 1 in the array holds sequence #0
% position 2 in the array holds sequence #1
% etc.
%
% Mike Lawrence 2009-04-03

%fprintf('PLEASE IGNORE THE FOLLOWING WARNINGS:\n'); 
%javaclasspath('/xchip/tcga/gbm/analysis/lawrence/samtools/trunk/sam-jdk/dist/sam-1.0.jar')
import net.sf.samtools.*;
import java.io.*;
import java.lang.*;
%fprintf('THANK YOU.\n\n');

r = SAMFileReader(File(bamfile));
s = r.getFileHeader.getSequences.iterator;
i = 1;
while(s.hasNext), D{i,1} = s.next.getSequenceName.toCharArray'; i=i+1; end

