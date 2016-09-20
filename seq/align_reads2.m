function [breakpoint,skip,score,seq2first,forseq,maxcounts]=align_reads2(dirnum,fastqname,seq1,seq2,seq_start,seq_end,startongenome,str1,str2)

enough_reads=4;
FASTQStruct = fastqread([dirnum '/' fastqname]);
[u,m]=unique({FASTQStruct.Header});
%if exist([dirnum '/' fastqname],'file')
    delete([dirnum '/' fastqname]);
%end
fastqwrite([dirnum '/' fastqname], FASTQStruct(m));
reads={FASTQStruct(m).Sequence};
isseq=cellfun(@length,reads)>30;
reads=reads(isseq);
n=length(reads);
revreads=cellfun(@seqrcomplement,reads,'UniformOutput',false);
breakpoints=[];
skips=[];
scores=[];
forseqs=[];
inversion=str1==str2;
seq1first=~str1;
for i=1:n
    l=length(reads{i});    
    if seq1first
        [score(1),match_only1(1),breakpoint1(1),breakpoint2(1),onread(1),skip(1)]=align_with_breakpoint(seq1, seq2, reads{i}, 4);
        [score(2),match_only1(2),breakpoint1(2),breakpoint2(2),onread(2),skip(2)]=align_with_breakpoint(seq1, seq2, revreads{i}, 4);
        [score ind]=max(score);
    else
        [score(1),match_only1(1),breakpoint2(1),breakpoint1(1),onread(1),skip(1)]=align_with_breakpoint(seq2, seq1, reads{i}, 4);
        [score(2),match_only1(2),breakpoint2(2),breakpoint1(2),onread(2),skip(2)]=align_with_breakpoint(seq2, seq1, revreads{i}, 4);
        [score ind]=max(score);
    end
    score=score/l;
    if ((~match_only1(ind))&&(score>0.8))
        breakpoint=([breakpoint1(ind);breakpoint2(ind)])+startongenome+seq_start-1+[0;skip(ind)];
        if inversion
            breakpoint(2)=startongenome(2)+seq_end(2)-breakpoint2(ind)+skip(ind);
        end
        breakpoints=[breakpoints; breakpoint'];
        skips=[skips;skip(ind)];
        scores=[scores;score];
        if skip(ind)>0
            if ind==2
                forseq=revreads{i}((onread(ind)-skip(ind)+1):onread(ind));
            else
                forseq=reads{i}((onread(ind)-skip(ind)+1):onread(ind));
            end
        else
            forseq='';
        end
        forseqs=[forseqs;{forseq}];
        if (length(scores)>=enough_reads)
            [tmp0,tmp0,n]=unique(breakpoints,'rows');
            counts=accumarray(n,1);
            if max(counts)>=enough_reads
                break;
            end
        end
    end
end
if (isempty(breakpoints))
    score = 0;
    skip = 0;
    breakpoint = 0;
    maxcounts = 0;
    seq2first = ~seq1first;
    forseq = '';
else
    [u,m,n]=unique([breakpoints skips],'rows');
    counts=accumarray(n,1);
    [s ind]=sort(counts,'descend');
    if (length(s)>=2)&&(s(1)==s(2))
        compares=sum(s==s(1));
        bkpntscr=zeros(1,compares);
        for j=1:compares
            bkpntscr(j)=mean(scores(n==ind(j)));
        end
        [m bk]=max(bkpntscr);
    else
        bk=1;
    end
    breakpoint=u(ind(bk),1:2);
    skip=u(ind(bk),3);
    score=mean(scores(n==ind(bk)));
    maxcounts=s(bk);
    seq2first=str1;
    pick1=find(n==ind(bk),1);
    forseq=forseqs{pick1};
end
