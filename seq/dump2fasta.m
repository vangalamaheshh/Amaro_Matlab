function dump2fasta(fasta,R,B,rnames)

trans='NNACGT';
for i=1:size(R,1)
    seq=B(R(i,9):(R(i,9)+R(i,5)-R(i,4)),1);
    miss=seq>60;
    seq(miss)=seq(miss)-64;   
    fprintf(fasta,'>%s\t%d\t%d\n%s\n',rnames(i).toCharArray',R(i,6),R(i,12),trans(seq(seq>-5)+2));
end

    