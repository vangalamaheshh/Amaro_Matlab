function [which,matchfreq,half]=best_swalign(seq,seqs,minscorediff,min2span)
% Yotam Drier, yotamd@gmail.com

l=length(seqs);
s=zeros(l,1);
as=cell(l,1);
ws=cell(l,1);
for i=1:l
    [s(i),as{i},ws{i}]=swalign(seq,seqs{i},'Alphabet','NT','GapOpen',20,'ExtendGap',20);
end
[score,which]=max(s);
half=-1;
rl=size(as{which},2);
matchfreq=mean(as{which}(2,:)=='|');
if (score-max(s(setdiff(1:l,which))) < minscorediff)||(rl<=min2span)
    which = -4;
    matchfreq = 0;
else
    %    if which > 6
    isrev=~mod(which,2);
    mid=length(seqs{which})/2;
    % min2span=10; %rl/4
    if (ws{which}(2)+rl>mid+min2span)&&(ws{which}(2)<mid-min2span)
        half=2; % read spans the (perhaps nonexisitng) breakpoint
    else
        half=~xor(ws{which}(2)<mid-min2span,isrev);
    end
    %end
end
