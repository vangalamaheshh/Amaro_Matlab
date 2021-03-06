addpath ~gadgetz/matlab
cd ~/projects/miRNA

c=read_colorscheme('~/projects/miRNA/data/colorscheme.txt');

cluster_and_display_res_cls('~/projects/miRNA/new/GCM/miGCM_218.gct',...
                            '~/projects/miRNA/new/GCM/new_cls_218.cls',...
                            struct('dist','correlation','preproc','row_center_and_normalize','cluster','average'),...
                            struct('dist','correlation','preproc','row_center_and_normalize','cluster','average'),...
                            'jun-reduced',c,'testout',{{'png','-r180'},{'pdf'}});

