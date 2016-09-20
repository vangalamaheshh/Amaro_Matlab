function N=view_allele_data(N, raw_bindirname, output_raw, S, Lsort, file_ext)
% 
% addpath ~/CancerGenomeAnalysis/trunk/matlab/
% addpath ~/CancerGenomeAnalysis/trunk/matlab/snp
% addpath ~/CancerGenomeAnalysis/trunk/matlab/POS
% addpath ~/CancerGenomeAnalysis/trunk/matlab/Startup

%make viewable bin directory for raw data
Nmin.dat=N.adat(:,:,1)
Nmin.chrn=N.chrn;
Nmin.chr=N.chr;
Nmin.chr=N.marker;
Nmin.chr=N.chr;
Nmin.marker=N.marker;
Nmin.pos=N.pos;
Nmin.sdesc=N.sdesc;
Nmin.sis=N.sis;
Nmax=Nmin;
Nmax.dat=N.adat(:,:,2);
% Ncom=Nmin;
% Ncom.dat=[Nmax.dat Nmin.dat];
% Ncomb=Ncom;
% for i=1:size(Nmin.dat,2) 
%     Ncom.sis(i).array=[Ncom.sis(i).array '_min'];
% end
% Ncomb.sis=[Ncom.sis N.sis];
% Ncomb.sdesc={Ncomb.sis.array}
if output_raw==1
write_bin([raw_bindirname '_min/'],Nmin)
write_bin([raw_bindirname '_max/'],Nmax)
end

%make viewable seg file for segmented data
Smin.dat=S.adat(:,:,1);
Smin.chrn=S.chrn;
Smin.chr=S.chr;
Smin.chr=S.marker;
Smin.chr=S.chr;
Smin.marker=S.marker;
Smin.pos=S.pos;
Smin.sdesc=S.sdesc;
Smin.sis=S.sis;
Smax=Smin;
Smax.dat=S.adat(:,:,2);
% Scomb=Smin;
% Scomb.dat=[Smax.dat Smin.dat]
% for i=1:size(Smin.dat,2) 
%     Scomb.sis(i).array=[Scomb.sis(i).array '_min']
% end
% Scomb.sis=[Scomb.sis S.sis]
% Scomb.sdesc={Scomb.sis.array};

write_seg_file([file_ext '.min.ascn.seg.txt'],Smin)
write_seg_file([file_ext '.max.ascn.seg.txt'],Smax)

%make viewable seg file for LOH data
write_seg_file([file_ext '.loh.seg.txt'],Lsort)