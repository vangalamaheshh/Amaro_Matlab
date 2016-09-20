function bsub_pipeline(sample)
ss = sample_to_short_sample(sample);

job1 =   bsub(['-R "rusage[matlab=1]" ''''matlab -nodisplay '...
         '-r "pipeline(''' sample ''')"'''''],[ss 'Pipeln']);
