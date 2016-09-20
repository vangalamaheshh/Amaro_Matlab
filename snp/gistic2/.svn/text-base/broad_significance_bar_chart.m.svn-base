function broad_significance_bar_chart(new_dir,broad_sig_file,D,min_qv,xlim_minqv)

cd(new_dir);
b=load_struct(broad_sig_file,'%s%f%f%f%f%f%f%f');

zA=b.Ampzscore;
zD=b.Delzscore;

disp('using one sided test');
%pA = 1-normcdf(zA,0,1);
pA = erfc(zA/sqrt(2))/2;
pA(pA<min_qv(1))=min_qv(1);

%pD = 1-normcdf(zD,0,1);
pD = erfc(zD/sqrt(2))/2;
pD(pD<min_qv(2))=min_qv(2);

qA = calc_fdr_value(pA);
qD = calc_fdr_value(pD);

qv_cutoff=0.25;
zA_range=[ max(zA(qA>=qv_cutoff)) min(zA(qA<qv_cutoff)) ];
zD_range=[ max(zD(qD>=qv_cutoff)) min(zD(qD<qv_cutoff)) ];

zA_cutoff=norminv(1-qv_cutoff/length(zA)*(length(find(qA<=qv_cutoff))+1),0,1);
zA_cutoff=clip_to_range(zA_cutoff,zA_range);

zD_cutoff=norminv(1-qv_cutoff/length(zD)*(length(find(qD<=qv_cutoff))+1),0,1);
zD_cutoff=clip_to_range(zD_cutoff,zD_range);

draw_arm_scores(D,b.Arm,qA,qv_cutoff,0,-1,[1 0 0],xlim_minqv(1));
outname='Amp.broad_significance_bar_chart';
saveas(gcf,[new_dir outname '.pdf'],'pdf');
saveas(gcf,[new_dir outname '_v3.eps'],'epsc2');
print_D([new_dir outname '_v2'],{{'pdf'},{'fig'}});


draw_arm_scores(D,b.Arm,qD,qv_cutoff,1,1,[0 0 1],xlim_minqv(2));
outname='Del.broad_significance_bar_chart';
saveas(gcf,[new_dir outname '.pdf'],'pdf');
saveas(gcf,[new_dir outname '_v3.eps'],'epsc2');
print_D([new_dir outname '_v2'],{{'pdf'},{'fig'}});



