function fuseseq3(bkpts_fn,startfrom,endat,expand_seq_around_bkpt,fastaname,refdir)
%
% Yotam Drier, yotamd@gmail.com
%
bkptsf=fopen(bkpts_fn);
bkpts=textscan(bkptsf,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s');
fclose(bkptsf);
if (~exist('endat','var') || endat == -1), endat=length(bkpts{1}); end
if ~exist('expand_seq_around_bkpt','var'), expand_seq_around_bkpt=100; end
if ~exist('fastaname','var'), fastaname='fused_seq3'; end
if ~exist('refdir','var'), refdir='hg18'; end
for idx=startfrom:endat
    b=bkpts{1}(idx);
    dirnum=num2str(b);
    if bkpts{9}(idx) > 0 % breakpoint found        
        if bkpts{10}(idx)==0 % if seq1 is first
            seq1=genome_region(bkpts{3}(idx),bkpts{4}(idx)-expand_seq_around_bkpt+1,bkpts{4}(idx),refdir);
            del1=genome_region(bkpts{3}(idx),bkpts{4}(idx)+1,bkpts{4}(idx)+expand_seq_around_bkpt,refdir);
            if bkpts{7}(idx) % inversion
                seq2=seqrcomplement(genome_region(bkpts{5}(idx),bkpts{6}(idx)-expand_seq_around_bkpt+1,bkpts{6}(idx),refdir));
                del2=seqrcomplement(genome_region(bkpts{5}(idx),bkpts{6}(idx)+1,bkpts{6}(idx)+expand_seq_around_bkpt,refdir));
            else
                seq2=genome_region(bkpts{5}(idx),bkpts{6}(idx)+1,bkpts{6}(idx)+expand_seq_around_bkpt,refdir);
                del2=genome_region(bkpts{5}(idx),bkpts{6}(idx)-expand_seq_around_bkpt+1,bkpts{6}(idx),refdir);
            end
        else
            if bkpts{7}(idx) % inversion
                seq1=seqrcomplement(genome_region(bkpts{5}(idx),bkpts{6}(idx)+1,bkpts{6}(idx)+expand_seq_around_bkpt,refdir));
                del1=seqrcomplement(genome_region(bkpts{5}(idx),bkpts{6}(idx)-expand_seq_around_bkpt+1,bkpts{6}(idx),refdir));
            else
                seq1=genome_region(bkpts{5}(idx),bkpts{6}(idx)-expand_seq_around_bkpt+1,bkpts{6}(idx),refdir);
                del1=genome_region(bkpts{5}(idx),bkpts{6}(idx)+1,bkpts{6}(idx)+expand_seq_around_bkpt,refdir);
            end
            seq2=genome_region(bkpts{3}(idx),bkpts{4}(idx)+1,bkpts{4}(idx)+expand_seq_around_bkpt,refdir);
            del2=genome_region(bkpts{3}(idx),bkpts{4}(idx)-expand_seq_around_bkpt+1,bkpts{4}(idx),refdir);
        end
        if bkpts{11}(idx)
            fusedseq=upper([seq1 char(bkpts{12}(idx)) seq2]);
        else
            fusedseq=upper([seq1 seq2]);
        end
        if bkpts{7}(idx) && bkpts{10}(idx)  % for competability
            fusedseq=seqrcomplement(fusedseq);
        end
        fid=fopen([dirnum '/' fastaname '.fasta'],'w');
        fprintf(fid,'>%d_%d_%d_%d_%d_%d\n%s\n',b,bkpts{2}(idx),bkpts{3}(idx),bkpts{4}(idx),bkpts{5}(idx),bkpts{6}(idx),fusedseq);
        fclose(fid);
        seq1=upper(seq1);
        del1=upper(del1);
        seq2=upper(seq2);
        del2=upper(del2);
        ls1=length(seq1);
        ls2=length(seq2);
        ld1=length(del1);
        ld2=length(del2);
        len=min(ls1,ld2);
        j1=0;
        while (j1<len)&&(seq1(end-j1)==del2(end-j1))
            j1=j1+1;
        end
        len=min(ls2,ld1);
        j2=0;
        while (j2<len)&&(del1(j2+1)==seq2(j2+1))
            j2=j2+1;
        end
        mh=j1+j2;                
        [score,a,w]=swalign([seq1(max(1,ls1-75):ls1) del1(1:min(ld1,75))],[del2(max(1,ld2-75):ld2) seq2(1:min(ls2,75))]);
        match=sum(a(2,:)=='|');
        if bkpts{11}(idx)==0 % when there's no insertion of non templated sequence, breakpoint is not well defined around the microhomoogy (if any), and so always define the largest (overlapping) segments
            if bkpts{10}(idx)==0 % if seq1 is first
                bp1=bkpts{4}(idx)+j2;
                if bkpts{7}(idx) % inversion
                    bp2=bkpts{6}(idx)+j1;
                else
                    bp2=bkpts{6}(idx)-j1;
                end
            else
                if bkpts{7}(idx) % inversion
                    bp2=bkpts{6}(idx)-j2;
                else
                    bp2=bkpts{6}(idx)+j2;
                end
                bp1=bkpts{4}(idx)-j1;
            end
        else
            bp1=bkpts{4}(idx);
            bp2=bkpts{6}(idx);
        end
        fid=fopen([dirnum '/' fastaname '.inf'],'w');
        fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n',b,bkpts{2}(idx),bkpts{3}(idx),bp1,bkpts{5}(idx),bp2,seq1,del1,bkpts{12}{idx},del2,seq2,mh,match,bkpts{10}(idx),bkpts{11}(idx));
        fclose(fid);
    else
        fid=fopen([dirnum '/' fastaname '.inf'],'w');
        fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n',b,bkpts{2}(idx),bkpts{3}(idx),bkpts{4}(idx),bkpts{5}(idx),bkpts{6}(idx),'failed','failed','failed','failed','failed',-1,-1,bkpts{10}(idx),bkpts{11}(idx));
        fclose(fid);
    end
end
