function [support,half]=classpair(seq,refs,minscorediff,align_min_qual,min2span)
% Yotam Drier, yotamd@gmail.com

support=-1;
half=-1;
if ~(strcmp(seq,'-1')||strcmp(seq,'[UNPAIRED_READ]'))
    l=length(refs);
    scores=zeros(l,1);
    alignments=cell(l,1);
    starts=cell(l,1);
    for i=1:l
        [scores(i),alignments{i},starts{i}]=swalign(seq,refs{i},'Alphabet','NT','GapOpen',20,'ExtendGap',20);
    end
    rls=cellfun('size',alignments,2);
    scores(rls<=2*min2span)=0;
    midfusion=length(refs{5})/2;
    if (starts{5}(2)+rls(5)<midfusion+min2span)||(starts{5}(2)>midfusion-min2span)
        scores(5)=0;
    end
    if (starts{6}(2)+rls(6)<midfusion+min2span)||(starts{6}(2)>midfusion-min2span)
        scores(6)=0;
    end
    [score,which]=max(scores);
    matchfreq=mean(alignments{which}(2,:)=='|');
    if (matchfreq > align_min_qual)&&(rls(which) > min2span)
        support=ceil(which/2);
        if (score-max(scores(setdiff(1:l,which))) < minscorediff)
            support = -2;
        else
            isrev=~mod(which,2);
            mid=length(refs{which})/2;
            if (starts{which}(2)+rls(which)>mid+min2span)&&(starts{which}(2)<mid-min2span)
                half=2; % read spans the (perhaps nonexisitng) breakpoint
            else
                half=~xor(starts{which}(2)<=mid-min2span,isrev);
                support = support+3;
            end
        end
    end
end
