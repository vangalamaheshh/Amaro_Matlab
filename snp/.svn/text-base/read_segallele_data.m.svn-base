function Smod=read_segallele_data(N, minoutputfile, maxoutputfile, totoutputfile, mindir, maxdir, totdir, havesegfile, n, method)
% Smod=read_segallele_data(N, minoutputfile, maxoutputfile, totoutputfile, mindir, maxdir, totdir, havesegfile, n, method)
% addpath ~/CancerGenomeAnalysis/trunk/matlab/
% addpath ~/CancerGenomeAnalysis/trunk/matlab/snp
% addpath ~/CancerGenomeAnalysis/trunk/matlab/POS
% addpath ~/CancerGenomeAnalysis/trunk/matlab/Startup

%  N=load_D2(infile)
%  N2=copyD(N);
%  N=convert_to_memfield(N, 'adat');

Min.dat=N.adat(:,:,1);
Min.chrn=N.chrn;
Min.chr=N.chr;
Min.pos=N.pos;
Max.dat=N.adat(:,:,2);
Max.chrn=N.chrn;
Max.chr=N.chr;
Max.pos=N.pos;
Max.sdesc=N.sdesc;
Min.sdesc=N.sdesc;

if havesegfile==0
    output=maxoutputfile;
    f=fopen(output,'w');
    for i=1:size(Max.dat,2)
        seg=read_dlm_file([ maxdir 'Sample' sprintf('%03d',i) '.1.01.seg.dat']);
        for j=1:length(seg)
            if seg{j}{3}~='N'
                seg{j}(2)=regexprep(seg{j}(2),'\"','');
                fprintf(f,['%d %s %s %s %s %s\n'],i,seg{j}{2},seg{j}{3},seg{j}{4},seg{j}{5},seg{j}{6});
                
            end
        end
    end
    fclose(f);
    output1=minoutputfile;
    f=fopen(output1,'w');
    for i=1:size(Min.dat,2)
        seg=read_dlm_file([ mindir 'Sample' sprintf('%03d',i) '.1.01.seg.dat']);
        for j=1:length(seg)
            if seg{j}{3}~='N'
                seg{j}(2)=regexprep(seg{j}(2),'\"','');
                fprintf(f,['%d %s %s %s %s %s\n'],i,seg{j}{2},seg{j}{3},seg{j}{4},seg{j}{5},seg{j}{6});
            end
        end
    end
    fclose(f);
    output2=totoutputfile;
    f=fopen(output2,'w');
    for i=1:size(N.dat,2)
        seg=read_dlm_file([ totdir 'Sample' sprintf('%03d',i) '.1.01.seg.dat']);
        for j=1:length(seg)
            if seg{j}{3}~='N'
                seg{j}(2)=regexprep(seg{j}(2),'\"','');
                fprintf(f,['%d %s %s %s %s %s\n'],i,seg{j}{2},seg{j}{3},seg{j}{4},seg{j}{5},seg{j}{6});
            end
        end
    end
    fclose(f);
    
    unix(['sed ''s/X/23/g'' ' minoutputfile ' > tempsegment_v1.txt']);
    unix([ 'sed ''s/Y/24/g'' tempsegment_v1.txt > tempsegment_v2.txt'])
    unix([ 'sort -n -k1 -k2 -k3 tempsegment_v2.txt > tempminsort.txt'])
    minoutputfile='tempminsort.txt';
    
    unix(['sed ''s/X/23/g'' ' maxoutputfile ' > tempsegment_v1.txt']);
    unix(['sed ''s/Y/24/g'' tempsegment_v1.txt > tempsegment_v2.txt'])
    unix(['sort -n -k1 -k2 -k3 tempsegment_v2.txt > tempmaxsort.txt'])
    maxoutputfile='tempmaxsort.txt';
    
    unix(['sed ''s/X/23/g'' ' totoutputfile ' > tempsegment_v1.txt']);
    unix(['sed ''s/Y/24/g'' tempsegment_v1.txt > tempsegment_v2.txt'])
    unix(['sort -n -k1 -k2 -k3 tempsegment_v2.txt > temptotsort.txt'])
    totoutputfile='temptotsort.txt';
end

Maxseg=read_cbs_file(maxoutputfile, Max, 0,[],0,1,0);
Minseg=read_cbs_file(minoutputfile, Min, 0, [],0,1,0);
Totseg=read_cbs_file(totoutputfile, N, 0, [],0,1,0);
set_verbose_level(10);

Minseg2=Minseg;
Maxseg2=Maxseg;
Totseg2=Totseg;
Minseg2=rmfield(Minseg2,'cbs_rl');
Maxseg2=rmfield(Maxseg2,'cbs_rl');
Totseg2=rmfield(Totseg2,'cbs_rl');

% Minseg2.cbs=2.^(Minseg.cbs+1);
% Maxseg2.cbs=2.^(Maxseg.cbs+1);
clear Min Max
%% modify segmentation
%keyboard
rl1=runlength(Minseg2.cbs,Minseg2.chrn); %min
rl2=runlength(Maxseg2.cbs,Maxseg2.chrn); %max

for i=1:length(rl1)
    start1{i}=rl1{i}(:,1);
    end1{i}=rl1{i}(:,2);
    start2{i}=rl2{i}(:,1);
    end2{i}=rl2{i}(:,2);
    combostart{i}=cat(1, start1{i}, start2{i});
    combostart{i}=sort(combostart{i});
    combostart{i}=unique(combostart{i});
    
    for j=1:length(combostart{i})-1
        if combostart{i}(j)==combostart{i}(j+1)
            if j+1==length(combostart{i})
                comboend{i}(j-1)=max(end1{i});
            else
                comboend{i}(j)=combostart{i}(j+2)-1;
            end
        else
            comboend{i}(j)=combostart{i}(j+1)-1;
        end
        comboend{i}(j+1)=max(end1{i});
    end
    
    for j=1:length(combostart{i})
        med1{i}(j)=nanmedian(N.adat(combostart{i}(j):comboend{i}(j),i,1));
        med2{i}(j)=nanmedian(N.adat(combostart{i}(j):comboend{i}(j),i,2));
    end
end

for i=1:length(rl1)
    rl11{i}(:,1)=combostart{i};
    rl11{i}(:,2)=comboend{i};
    rl11{i}(:,3)=med1{i};
    rl21{i}(:,1)=combostart{i};
    rl21{i}(:,2)=comboend{i};
    rl21{i}(:,3)=med2{i};
end

Minseg2.cbs=derunlength(rl11);
Maxseg2.cbs=derunlength(rl21);
clear rl1 rl2  start1 start2 end1 end2  combostart comboend med1 med2
clear Minseg Maxseg Totseg
%need to filter combined list of breakpoints against total cn breakpoints
for i=1:length(rl11)
    for j=1:length(rl11{i})-1
        ind1{i}=find((rl11{i}(j,3)+rl21{i}(j,3))==(rl11{i}(j+1,3)+rl21{i}(j+1,3))+.1 );
        ind2{i}=find((rl11{i}(j,3)+rl21{i}(j,3))==(rl11{i}(j+1,3)+rl21{i}(j+1,3))-.1);
        if ~isempty(ind1{i} | ind2{i})
            warning('allele-specific copy number changes from 2,0 to 1,1 (or reverse) - should not smooth over segment')
        end
    end
end

S.sdesc=N.sdesc;
S.sis=N.sis;

rl3=runlength(Totseg2.cbs,Totseg2.chrn);

switch method
    case 'combine'
        for i=1:length(rl11)
            a{i}(:,1)=rl11{i}(:,1);
            a{i}(:,2)=rl11{i}(:,2);
            a{i}(:,3)=1; %allele
            t{i}(:,1)=rl3{i}(:,1);
            t{i}(:,2)=rl3{i}(:,2);
            t{i}(:,3)=0; %total
            totallele{i}=cat(1, a{i}, t{i});
            [totallele{i}(:,1),IX]=sort(totallele{i}(:,1));
            totallele{i}(:,2)=totallele{i}(IX,2);
            totallele{i}(:,3)=totallele{i}(IX,3);
            [b,m,m1]=unique(totallele{i}(:,1),'rows');
            ta{i}(:,1)=totallele{i}(m,1);
            ta{i}(:,2)=totallele{i}(m,2);
            ta{i}(:,3)=totallele{i}(m,3);
            
            for j=1:length(ta{i})-1
                if ta{i}(j,1)==ta{i}(j+1,1)
                    if j+1==length(ta{i})
                        ta{i}(j-1,2)=max(a{i}(:,2));
                    else
                        ta{i}(j,2)=ta{i}(j+2,1)-1;
                    end
                else
                    ta{i}(j,2)=ta{i}(j+1,1)-1;
                end
                ta{i}(j+1,2)=max(a{i}(:,2));
            end
            % indtot{i}=find(ta{i}(:,3)==0);
            % indallele{i}=find(ta{i}(:,3)==1);
        end
        clear m m1 t a b IX
        Minseg2.dat=derunlength(ta);
        rl4=runlength(Minseg2.dat,Minseg2.chrn);
        for i=1:length(rl4)
            for j=1:length(rl4{i})
                rl1a{i}(j,1)=rl4{i}(j,1);
                rl1a{i}(j,2)=rl4{i}(j,2);
                rl1a{i}(j,3)=nanmedian(N.adat(rl4{i}(j,1):rl4{i}(j,2),i,1));
                rl2a{i}(j,1)=rl4{i}(j,1);
                rl2a{i}(j,2)=rl4{i}(j,2);
                rl2a{i}(j,3)=nanmedian(N.adat(rl4{i}(j,1):rl4{i}(j,2),i,2));
            end
        end
        Minseg2.cbs=derunlength(rl1a);
        Maxseg2.cbs=derunlength(rl2a);
        
    case 'filter' %only total cn segments
        for i=1:length(rl3)
            for j=1:length(rl3{i})
                rl1a{i}(j,1)=rl3{i}(j,1);
                rl1a{i}(j,2)=rl3{i}(j,2);
                rl1a{i}(j,3)=nanmedian(N.adat(rl3{i}(j,1):rl3{i}(j,2),i,1));
                rl2a{i}(j,1)=rl3{i}(j,1);
                rl2a{i}(j,2)=rl3{i}(j,2);
                rl2a{i}(j,3)=nanmedian(N.adat(rl3{i}(j,1):rl3{i}(j,2),i,2));
            end
        end
        Minseg2.cbs=derunlength(rl1a);
        Maxseg2.cbs=derunlength(rl2a);
end

clear rl11 rl21 rl1a rl2a rl3

%% smooth
%n= number of SNPs to smooth
if max(Minseg2.chrn)==23 || max(Minseg2.chrn)==24
[noXY,x1,x2]=match_string_sets({'23'},Minseg2.chr);
else
[noXY,x1,x2]=match_string_sets({'X'},Minseg2.chr);
end

if ~isempty(x2)
Minseg2_noXY=reorder_D_rows(Minseg2,1:(x2-1));
Maxseg2_noXY=reorder_D_rows(Maxseg2,1:(x2-1));
Totseg2_noXY=reorder_D_rows(Totseg2,1:(x2-1));
else
    Minseg2_noXY=Minseg2;
Maxseg2_noXY=Maxseg2;
Totseg2_noXY=Totseg2;
end

clear Totseg2 Minseg2 Maxseg2

Totseg2_noXY.cbs=(2.^(Totseg2_noXY.cbs)+1);

 
Totseg3=smooth_cbs2(Totseg2_noXY,n);
Minseg3=smooth_cbs2(Minseg2_noXY,n);
Maxseg3=smooth_cbs2(Maxseg2_noXY,n);

Smod=S;
Smod.chr=Totseg3.chr;
Smod.pos=Totseg3.pos;
Smod.chrn=Totseg3.chrn;
Smod.marker=Totseg3.marker;
Smod.adat(:,:,1)=Minseg3.cbs;
Smod.adat(:,:,2)=Maxseg3.cbs;
Smod.dat=Totseg3.cbs;

%% convert Smod to hdf5 and save for use in subsequent functions
% DS=datastruct(Smod);
% DS=convert_to_diskfield(DS, 'adat', 'allelic.segmented.dat.h5')
% save_D2(outfile,DS);
