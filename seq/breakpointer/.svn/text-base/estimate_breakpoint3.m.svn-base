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
%   omitted 11. the length of micro homology
%   11. the length of foreign sequence added at the breakpoint
%   12. the foreign sequence added at the breakpoint

function estimate_breakpoint3(indlist,fileid,fastqname,params)

if ~exist('params','var'), params=[]; end
if ~isfield(params,'fish_no_supp_reads_thres'), params.fish_no_supp_reads_thres=6; end
if ~isfield(params,'fish_low_confidence_fixsidewithread'), params.fish_low_confidence_fixsidewithread=80; end
if ~isfield(params,'fish_low_confidence_fixsidewithoutread'), params.fish_low_confidence_fixsidewithoutread=200; end
if ~isfield(params,'fish_high_confidence_fixsidewithread'), params.fish_high_confidence_fixsidewithread=40; end
if ~isfield(params,'fish_high_confidence_fixsidewithoutread'), params.fish_high_confidence_fixsidewithoutread=50; end
if ~isfield(params,'align_enough_reads'), params.align_enough_reads=20; end
if ~isfield(params,'align_wsplit'), params.align_wsplit=8; end
if ~isfield(params,'align_min_qual'), params.align_min_qual=0.75; end
if ~isfield(params,'max_read_length'), params.max_read_length=150; end
outputf=fopen(['my_breakpoints_' fileid '.txt'],'w');
%warning('OFF','Bioinfo:swalign:EmptyAlignment');
for idx=1:length(indlist)
    b=indlist(idx);
    if mod(idx,10)==1
        [idx b]
    end
    dirnum=num2str(b);
    fid=fopen([dirnum '/splitreads.helper2']);
    t=fscanf(fid,'%d\t',11);
    seq1=fscanf(fid,'%s\n',1);
    seq2=fscanf(fid,'%s\n',1);
    fclose(fid);
    inversion=t(4)==t(7);
    seq1=upper(seq1);
    seq2=upper(seq2);
    chr=[t(2);t(5)];
    startongenome=[t(3);t(6)]-t(9);
    if t(8)<=params.fish_no_supp_reads_thres
        fixsidewithread=params.fish_low_confidence_fixsidewithread;
        fixsidewithoutread=params.fish_low_confidence_fixsidewithoutread;
    else
        fixsidewithread=params.fish_high_confidence_fixsidewithread;
        fixsidewithoutread=params.fish_high_confidence_fixsidewithoutread;
    end
    fix=[fixsidewithread fixsidewithoutread];
    breakpointstart=[t(3);t(6)]-[fix(t(4)+1);fix(t(7)+1)];
    breakpointend=[t(3);t(6)]+[fix(2-t(4));fix(2-t(7))];
    bprlen=breakpointend-breakpointstart+1;
    bpoffset=breakpointstart-startongenome;
    seq_start=max(bpoffset-params.max_read_length,1);
    seq_end=min(bpoffset+bprlen+params.max_read_length,([length(seq1); length(seq2)]));
    seq1=seq1(seq_start(1):seq_end(1));
    seq2=seq2(seq_start(2):seq_end(2));
    if inversion
        seq2=seqrcomplement(seq2);
    end
    [breakpoint,skip,score,seq2first,forseq,maxcounts]=align_reads3(dirnum,fastqname,seq1,seq2,seq_start,seq_end,startongenome,t(4),t(7),params);
    if (score>0)
        fprintf(outputf,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n',b,t(1),chr(1),breakpoint(1),chr(2),breakpoint(2),inversion,maxcounts,score,seq2first,skip,forseq);
    else
        fprintf(outputf,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n',b,t(1),chr(1),-1,chr(2),-1,-1,-1,-1,-1,-1,'failed');
    end
end
fclose(outputf);
%warning('ON','Bioinfo:swalign:EmptyAlignment');

