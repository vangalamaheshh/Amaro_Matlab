function dump2fastq(fastq,R,B,rnames)

trans='NNACGT';
for i=1:size(R,1)
    seq=B(R(i,9):(R(i,9)+R(i,5)-R(i,4)),1);
    bc=B(R(i,9):(R(i,9)+R(i,5)-R(i,4)),2);
    miss=seq>60;
    seq(miss)=seq(miss)-64;  
    fprintf(fastq,'@%s\n%s\n+\n%s\n',rnames(i).toCharArray',trans(seq(seq>-5)+2),char(bc(bc>-5)+64));
end

    
