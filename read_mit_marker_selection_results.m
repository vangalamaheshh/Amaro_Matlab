function [r,fet,s,p,fpr,fwer,rankp,fdr,q,pi0,T]=read_mit_marker_selection_results(fname)

T=read_mit_odf_file(fname);
r=T.data{1};
fet=strvcat(T.data{2});
s=T.data{3};
p=T.data{4};
fpr=T.data{5};
fwer=T.data{6};
rankp=T.data{7};
fdr=T.data{8};
q=T.data{9};
pi0=str2num(T.pi0);
