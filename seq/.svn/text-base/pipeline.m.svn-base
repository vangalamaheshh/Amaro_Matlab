function pipeline(sample)

if iscell(sample), error('multiple samples not supported.'); end

bd1 = '/xchip/cga1/firehose_output/Individual';
bd2 = '/xchip/tcga_scratch/lawrence';
bd3 = '/xchip/tcga_scratch/ng';

ss = sample_to_short_sample(sample);
name1 = upper(regexprep(sample,'/','-'));
name1 = regexprep(name1,'-WGS$','/wgs');
name1 = regexprep(name1,'-CAPTURE$','/capture');
name2 = lower(regexprep(sample,'-','/'));

dir1 = [bd1 '/' name1];
dir2 = [bd2 '/' name2];
dir3 = [bd3 '/' name1];

if ~exist(dir2,'dir'), error('Not found: %s',dir2); end

job1 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "make_lanetable(''' sample ''')"'''''],[ss 'LaneT']);
job2 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "make_lanelist(''' sample ''')"'''''],[ss 'LaneL']);
job3 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "CoverageByBase(''' sample ''')"'''''],[ss 'CovBB']);
job4 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "InsertSizeByLane(''' sample ''')"'''''],[ss 'ISize']);
job5 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "SegSeqPreprocess(''' sample ''')"'''''],[ss 'SSPre']);
job6 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "dRangerPreprocess(''' sample ''')"'''''],[ss 'dRPre']);

% coverage stats
bwait(job3);
job7 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "get_global_coverage_stats(''' sample ''',''global_coverage_by_zone.txt'',''zone'')"'''''],[ss 'Zone']);
job8 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "extract_from_cbb(''' sample ''',' ...
   '''/xchip/tcga_scratch/lawrence/capture/whole_exome_refseq_coding.targets.interval_list.GENESGC.txt'','...
   '''we_genes_category_coverage.txt'','...
   '''/xchip/tcga_scratch/lawrence/db/context'',4)"'''''],[ss 'Context']);

% SegSeq
bwait(job5);

dd = [dir3 '/cn'];
if ~exist(dd,'dir')
    fprintf('Creating: %s\n',dd);
    mkdir(dd);
end
system(['ln -s ' dir2 '/tumor_SS ' dd ';'...
        'ln -s ' dir2 '/normal_SS ' dd]);
job9 =  bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "SegSeqRun(''' sample ''')"'''''],[ss 'SegSeq']);

% dRanger
bwait([job1;job2;job6;job4]);
job10 =   bsub(['-q hugemem -R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "dRangerRun(''' sample ''')"'''''],[ss 'dRanger']);
