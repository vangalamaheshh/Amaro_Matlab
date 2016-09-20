function fh_BreakPointerScatterTest(sample,drr,bam,blacklist,refdir,insertionsize_input,fish_no_supp_reads_thres,fish_low_confidence_fixsidewithread,fish_low_confidence_fixsidewithoutread,fish_high_confidence_fixsidewithread,fish_high_confidence_fixsidewithoutread,fish_allowedmm,fish_tipsize,fish_minmmintip,fish_maxN,fish_expand_pairs_extraction,fish_maxreads,readlen,align_enough_reads,align_wsplit,align_min_qual,libdir,firstidx,lastidx)
% Mike Lawrence 2010-02-05
%
% wraps BreakPointer by Yotam Drier 2010

if nargin~=24, error('Usage: fh_BreakPointerScatter(sample,drr,bam,blacklist,refdir,insertionsize_input,fish_no_supp_reads_thres,fish_low_confidence_fixsidewithread,fish_low_confidence_fixsidewithoutread,fish_high_confidence_fixsidewithread,fish_high_confidence_fixsidewithoutread,fish_allowedmm,fish_tipsize,fish_minmmintip,fish_maxN,fish_expand_pairs_extraction,fish_maxreads,readlen,align_enough_reads,align_wsplit,align_min_qual,libdir,firstidx,lastidx)'); end

if strcmpi(sample,'NO_DATA')
  fprintf('fh_BreakPointerScatter exiting because received "NO_DATA" signal\n');
  return;
end

demand_file(drr); demand_file(bam);
if ~strcmp(blacklist,'none'), demand_file(blacklist); end
if ~isnumeric(fish_no_supp_reads_thres), fish_no_supp_reads_thres=str2double(fish_no_supp_reads_thres); end
if ~isnumeric(fish_allowedmm), fish_allowedmm=str2double(fish_allowedmm); end
if ~isnumeric(fish_tipsize), fish_tipsize=str2double(fish_tipsize); end
if ~isnumeric(fish_minmmintip), fish_minmmintip=str2double(fish_minmmintip); end
if ~isnumeric(fish_maxN), fish_maxN=str2double(fish_maxN); end
if ~isnumeric(fish_low_confidence_fixsidewithread), fish_low_confidence_fixsidewithread=str2double(fish_low_confidence_fixsidewithread); end
if ~isnumeric(fish_low_confidence_fixsidewithoutread), fish_low_confidence_fixsidewithoutread=str2double(fish_low_confidence_fixsidewithoutread); end
if ~isnumeric(fish_high_confidence_fixsidewithread), fish_high_confidence_fixsidewithread=str2double(fish_high_confidence_fixsidewithread); end
if ~isnumeric(fish_high_confidence_fixsidewithoutread), fish_high_confidence_fixsidewithoutread=str2double(fish_high_confidence_fixsidewithoutread); end
if ~isnumeric(align_enough_reads), align_enough_reads=str2double(align_enough_reads); end
if ~isnumeric(align_wsplit), align_wsplit=str2double(align_wsplit); end
if ~isnumeric(align_min_qual), align_min_qual=str2double(align_min_qual); end
if ~isnumeric(readlen), readlen=str2double(readlen); end
if ~isnumeric(fish_expand_pairs_extraction), fish_expand_pairs_extraction=str2double(fish_expand_pairs_extraction); end
if ~isnumeric(fish_maxreads), fish_maxreads=str2double(fish_maxreads); end
if ~isnumeric(insertionsize_input)
    if isempty(strfind(insertionsize_input,'.isz'))
        if isempty(strfind(insertionsize_input,','))
            maxinsertionsize = str2double(insertionsize_input);
            mininsertionsize = maxinsertionsize;
        else
            is=str2double(regexp(insertionsize_input,',','split'));
            mininsertionsize=is(1);
            maxinsertionsize=is(2);
        end
    else
        demand_file(insertionsize_input);
        inshist = dlmread(insertionsize_input, '\t');
        %[m i]=max(sum(inshist(:,2:end),2));
        t=sum(inshist(:,2:end),2);
        [s,i]=sort(t,'descend');
        a=cumsum(s);
        third=find(a>a(end)*.75,1,'first');
        maxinsertionsize=max(i(1:third));
        mininsertionsize=min(i(1:third));
    end
else
    mininsertionsize=insertionsize_input;
    maxinsertionsize=insertionsize_input;
end
if ~isdir(refdir), error('Not found: reference directory %s',refdir); end
if ~isnumeric(firstidx), firstidx=str2double(firstidx); end
if ~isnumeric(lastidx), lastidx=str2double(lastidx); end
if firstidx<1 || lastidx<1 || firstidx>lastidx, error('problem with firstidx,lastidx'); end
aeb_script = fullfile(libdir,'align_each_bkpt5.pl'); demand_file(aeb_script);
gsr_jar = fullfile(libdir,'GrabSplitReads.jar'); demand_file(gsr_jar);

fileid = 'all';
fileid2 = 'bp';
fusedname = 'fused_seq';
splitname = 'splitreads.fastq';

fprintf('prepare temp files\n');
P=[];
P.refdir = refdir;
P.blacklist = blacklist;
P.gsr_jar = gsr_jar;
P.fastq_name = splitname;
P.maxreads = fish_maxreads;
P.fish_allowedmm=fish_allowedmm;
P.fish_tipsize=fish_tipsize;
P.fish_minmmintip=fish_minmmintip;
P.fish_maxN=fish_maxN;
P.readlen=readlen;
P.expand_pairs_extraction=fish_expand_pairs_extraction;
P.expand_seq_around_bkpt=fish_expand_pairs_extraction;
P.max_insertion_size = maxinsertionsize;
if exist('mininsertionsize','var'), P.min_insertion_size = mininsertionsize; end

prepare_breakpointer_temp_files3XXX(drr,bam,firstidx,lastidx,P);

fprintf('estimate_breakpoint\n');
params=struct('fish_no_supp_reads_thres',fish_no_supp_reads_thres,'fish_low_confidence_fixsidewithread',fish_low_confidence_fixsidewithread,'fish_low_confidence_fixsidewithoutread',fish_low_confidence_fixsidewithoutread,'fish_high_confidence_fixsidewithread',fish_high_confidence_fixsidewithread,'fish_high_confidence_fixsidewithoutread',fish_high_confidence_fixsidewithoutread,'align_enough_reads',align_enough_reads,'align_wsplit',align_wsplit,'align_min_qual',align_min_qual);
estimate_breakpoint3(firstidx:lastidx,fileid,splitname,params);

fprintf('fuseseq\n');
fname = ['my_breakpoints_' fileid '.txt'];
fuseseq3(fname,1,-1,100,fusedname,refdir);

fprintf('align_each_bkpt\n');
system(['perl ' aeb_script ' ' num2str(firstidx) ' ' num2str(lastidx) ' ' fusedname '.fasta ' fileid2 ' ' splitname ' ' libdir]);

fprintf('done\n');
