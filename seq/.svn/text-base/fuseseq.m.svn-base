function fuseseq(bkpts_fn,startfrom,endat)

expand_seq_around_bkpt=100;
bkptsf=fopen(bkpts_fn);
bkpts=textscan(bkptsf,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s');    
fclose(bkptsf);
if ~exist('endat','var'), endat=length(bkpts{1}); end
for idx=startfrom:endat
    b=bkpts{1}(idx);

    dirnum=num2str(b);
    fid=fopen([dirnum '/splitreads.helper']);
    t=fscanf(fid,'%d\t',9);
    seq1=fscanf(fid,'%s\n',1);
    seq2=fscanf(fid,'%s\n',1);
    fclose(fid);   
    bkptno=t(1);
    startongenome=[t(3);t(6)]-t(9);
    breakpointstart=[t(3);t(6)]-40;
    breakpointend=[t(3);t(6)]+50;
    if t(8)<=9
        breakpointstart=breakpointstart-40;
        breakpointend=breakpointend+340;
    end
    bprlen=breakpointend-breakpointstart+1;
    bpoffset=breakpointstart-startongenome;

    inversion=t(4)==t(7);
    breakpoint=[bkpts{4}(idx);bkpts{6}(idx)]-startongenome+1-[0;bkpts{12}(idx)]; %-seq_start
    if inversion
        seq2=seqrcomplement(seq2);
        breakpoint(2)=startongenome(2)-bkpts{6}(idx)+bkpts{12}(idx)+length(seq2)+1; %+seq_end(2)
    end       
    if bkpts{10}(idx)
        fusedseq=[seq2(max(breakpoint(2)-expand_seq_around_bkpt,1):breakpoint(2)) char(bkpts{13}(idx)) seq1(breakpoint(1):min(breakpoint(1)+expand_seq_around_bkpt,length(seq1)))];
    else
        fusedseq=[seq1(max(breakpoint(1)-expand_seq_around_bkpt,1):breakpoint(1)) char(bkpts{13}(idx)) seq2(breakpoint(2):min(breakpoint(2)+expand_seq_around_bkpt,length(seq2)))];
    end
    if (inversion)&&(bkpts{10}(idx))
        fusedseq=seqrcomplement(fusedseq);
    end

    fid=fopen([dirnum '/fused_seq.fasta'],'w');
    fprintf(fid,'>%d_%d_%d_%d_%d_%d\n%s\n',b,bkptno,bkpts{3}(idx),bkpts{4}(idx),bkpts{5}(idx),bkpts{6}(idx),fusedseq);
    fclose(fid);
end
