function display_figure(matname,fname)

load(matname);
display_D(D1ord,gdend,sdend,'jun-reduced');
print_D(fname,{{'png','-r180'}},1,[]);

