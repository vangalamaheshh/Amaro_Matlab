function [CNsort, Lsort, D]=LOH_byallele(S, mastersampleinfo, ignore_useme, arraylist)

% S=load_D2(infile)
% S2=copyD(S);
% S=convert_to_memfield(S, 'adat');
% 
% N=load_D2(infile2)
% N2=copyD(N);
% N=convert_to_memfield(N, 'adat');

%LOH calls on segmented data
S.loh(1:length(S.marker),1:size(S.adat,2))=NaN;
M.loh=((S.adat(:,:,1)<0.5 & S.adat(:,:,2)>1.5) & S.adat(:,:,2)<2.3);
S.loh(M.loh)=3; %copy neutral
M.loh=((S.adat(:,:,1)<0.5 & S.adat(:,:,2)<1.5));
S.loh(M.loh)=1; %copy loss
M.loh=((S.adat(:,:,1)>0.5 & S.adat(:,:,2)>0.5) & S.adat(:,:,2)<1.5);
S.loh(M.loh)=2; %retention
M.loh=((S.adat(:,:,1)<0.5 & S.adat(:,:,2)>2.3));
S.loh(M.loh)=4; %amp with LOH
M.loh=((S.adat(:,:,1)>0.8 & S.adat(:,:,2)>1.5))
S.loh(M.loh)=5; %amp with retention

D.dat=S.loh;
D.marker=S.marker;
D.pos=S.pos;
D.chr=S.chr;
D.chrn=S.chrn;
D.sis=S.sis
D.sdesc={S.sis.array};

MSI=read_sample_info_file(mastersampleinfo);

% remove normals
type1=strmatch({'T'},{D.sis.tumornormal});

% match 'use me'
if ignore_useme==0;
[U, u1, u2]=match_string_sets_hash({MSI.array},{D.sis.array});
use=MSI(u1);
type2=strmatch({'x'},{use.useme});
else
    type2=type1;
end
type=intersect(type1,type2);

if exist('arraylist', 'var')
    AL=read_array_list_file(arraylist)
    use_arrays={AL.array};
    [UD,idx,d2]=intersect({D.sis.array},use_arrays);
    if length(UD)~=length(use_arrays)
        warning('Did not match all arrays')
    end
    filter=intersect(type,idx);
else
    filter=type;
end

% D.sis=D.sdesc;
% D.sdesc={D.sis.array}
%% copy neutral LOH calls for gistic
L=D;
M.loh=((S.adat(:,:,1)<0.5 & S.adat(:,:,2)>1.5) & S.adat(:,:,2)<2.3);
L.dat(1:length(S.marker),1:size(S.loh, 2))=0;
L.dat(M.loh)=1;
% filter=intersect(type, good_samples);
CNsort=reorder_D_cols(L,filter);
clear M L
%% LOH calls for gistic
L1=D;

M.loh=((S.adat(:,:,1)<0.5 & S.adat(:,:,2)>1.5) & S.adat(:,:,2)<2.3);
M1.loh=((S.adat(:,:,1)<0.5 & S.adat(:,:,2)<1.5));
L1.dat(1:length(S.marker),1:size(S.loh, 2))=0;
L1.dat(M.loh)=1;
L1.dat(M1.loh)=1;

% filter=intersect(type, good_samples);
Lsort=reorder_D_cols(L1,filter);
clear M1 M L1
%% plots
display_D(Lsort,[],[],'lohsup');
print_D('all_lohcalls',{{'pdf'},{'png','-r180'}});

Dsort=reorder_D_cols(D,filter);
display_D(Dsort,[],[],'snpsup');
h=get(gcf,'Children');
ampmap=[0.8,0.8,0.8; 0.1,0.3,0.5; 1,1,0.4; 0,0.5,0.5; 1,0.6,0.2; 1,0.2,0.2]
subplot(h(4));
colormap(ampmap)
subplot(h(4));
caxis([0 5.5])
% colorbar;
print_D('allelecalls',{{'pdf'},{'png','-r180'}});

display_D(CNsort,[],[],'lohsup');
print_D('cn_lohcalls',{{'pdf'},{'png','-r180'}});

%% convert CNsort and Lsort to hdf5 and save for use in subsequent functions
% CN=datastruct(CNsort);
% CN=convert_to_diskfield(CNsort, 'dat', 'copyneutral.LOH.dat.h5')
% save_D2(outfileCN,CN);
% 
% L=datastruct(Lsort);
% L=convert_to_diskfield(Lsort, 'dat', 'all.LOH.dat.h5')
% save_D2(outfileLOH,L);