function res=test_diff_anal_step(D,cls0,cls1)

res=cell(2,1);
nparts=1;
verbose('1')
[P,S]=differential_analysis(D,cls0,cls1,'ttest');
verbose('2')
[PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1, ...
                                                  'ttest',5000,0,nparts,'/xchip/data/gadgetz/lsfres/');
verbose('3')
[Pf,PRf]=fix_Pvalues_by_parts(P,PR,nparts,'/xchip/data/gadgetz/lsfres/');

res{1}=Pf;
verbose('4')

[PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1, ...
                                                  'ttest',5000,1,nparts,'/xchip/data/gadgetz/lsfres/');
verbose('5')
[Pf,PRf]=fix_Pvalues_by_parts(P,PR,nparts,'/xchip/data/gadgetz/lsfres/');
verbose('6')

res{2}=Pf;


