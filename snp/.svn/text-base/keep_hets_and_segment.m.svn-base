function [N1,G1]=keep_hets_and_segment(N,G, seg, mindir, maxdir, totdir, arraylist)
% [N1,G1]=keep_hets_and_segment(N,G, seg, mindir, maxdir, totdir, arraylist)
% added segmentation to this function also, 
%if seg=1 - segment, if seg=0 - don't segment
% addpath ~/CancerGenomeAnalysis/trunk/matlab/
% addpath ~/CancerGenomeAnalysis/trunk/matlab/snp
% addpath ~/CancerGenomeAnalysis/trunk/matlab/POS
% addpath ~/CancerGenomeAnalysis/trunk/matlab/Startup

if isempty(seg)
    seg=1;
    if isempty(mindir)
        mindir='./'
    end
    if isempty(maxdir)
        maxdir='./'
    end
end

%N and G are both ordered by position
M=N;
%pull out hets for each tumor based on paired normal, NaNs for homozygous SNPs
if isfield(N.sis, 'collaboratorparticipantid')
for i=1:length(M.sdesc)
    if strcmp({N.sis(i).tumornormal},'Normal')==1
        p{i}=strmatch({N.sis(i).collaboratorparticipantid},{N.sis.collaboratorparticipantid});
        [snpr,scol,v]=find(G.dat(:,i) ~= 1);
        M.adat(snpr,p{i},:)=NaN;
    else p{i}=0;
    end
end
else %for old sample info with 'paired' column
 for i=1:length(M.sdesc)
    if strcmp({N.sis(i).tumornormal},'Normal')==1
        p{i}=strmatch({N.sis(i).name},{N.sis.paired});
        [snpr,scol,v]=find(G.dat(:,i) ~= 1);
        M.adat(snpr,p{i},:)=NaN;
    else p{i}=0;
    end
 end
end


%optional arraylist
if exist('arraylist', 'var')
    AL=read_array_list_file(arraylist)
    use_arrays={AL.array};
    [UD,idx,d2]=intersect({M.sis.array},use_arrays);
    if length(UD)~=length(use_arrays)
        warning('Did not match all arrays')
    end
    M=reorder_D_cols(M, idx);
  
    G1=reorder_D_cols(G,idx);
else
    G1=G;
end


N.adat=M.adat;
N.marker=M.marker;
N.pos=M.pos;
N.chr=M.chr;
N.chrn=M.chrn;
N.sdesc=M.sdesc;
N.sis=M.sis;
N.dat=M.dat;

% segment with cbs
if seg==1
    Min.dat=N.adat(:,:,1);
    Min.chrn=N.chrn;
    Min.chr=N.chr;
    Min.pos=N.pos;
    mkdir(mindir)
    cd(mindir)
    run_cbs_batch(Min,[],1,0,100);
    %log and don't smooth
    %unix('rm *1.dat')

    Max.dat=N.adat(:,:,2);
    Max.chrn=N.chrn;
    Max.chr=N.chr;
    Max.pos=N.pos
    mkdir(maxdir)
    cd(maxdir)
    run_cbs_batch(Max,[],1,0,100);
    %unix('rm *1.dat')

    mkdir(totdir)
    cd(totdir)
    run_cbs_batch(N,[],1,0,100);
    %unix('rm *1.dat')
end
cd ..

N=rmfield(N,'min');
N=rmfield(N,'max');
N1=N;

