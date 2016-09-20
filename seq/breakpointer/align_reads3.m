function [breakpoint,skip,score,seq2first,forseq,maxcounts]=align_reads3(dirnum,fastqname,seq1,seq2,seq_start,seq_end,startongenome,str1,str2,params)
%
% Yotam Drier, yotamd@gmail.com
%
breakpoint=0;
skip=0;
score=0;
maxcounts=0;
seq2first=0;
forseq='';
if exist([dirnum '/' fastqname],'file')
    FASTQStruct = fastqread([dirnum '/' fastqname]);
    reads={FASTQStruct.Sequence};
    %qual=cellfun(@(x) x-33,{FASTQStruct(m).Quality},'UniformOutput',false);
    qual={FASTQStruct.Quality};
    isseq=cellfun(@length,reads)>30;
    reads=reads(isseq);
    qual=qual(isseq);
    n=length(reads);
    revreads=cellfun(@seqrcomplement,reads,'UniformOutput',false);
    breakpoints=[];
    skips=[];
    scores=[];
    forseqs=[];
    inversion=str1==str2;
    seq1first=~str1;
    ls1=length(seq1);
    ls2=length(seq2);
    for i=1:n
        p=1-10.^(-(qual{i}-33)/10);
        if seq1first
            [score(1),match_only1(1),breakpoint1(1),breakpoint2(1),onread(1),skip(1)]=align_with_breakpoint(seq1, seq2, reads{i}, qual{i}, params.align_wsplit);
            % reads{i}(1:onread(1)-skip(1))
            % seq1(breakpoint1(1)-onread(1)+1+skip(1):breakpoint1(1))
            % reads{i}(onread(1)+1-skip(1):end)
            % seq2(breakpoint2(1)+1:breakpoint2(1)+length(reads{i})-onread(1))
            [score(2),match_only1(2),breakpoint1(2),breakpoint2(2),onread(2),skip(2)]=align_with_breakpoint(seq1, seq2, revreads{i}, qual{i}(end:-1:1), params.align_wsplit);
            % revreads{i}(1:onread(2)-skip(2))
            % seq1(breakpoint1(2)-onread(2)+1+skip(2):breakpoint1(2))
            % seq2(breakpoint2(2)+1:breakpoint2(2)+length(reads{i})-onread(2))==reads{i}(onread(2)+1:end)
            [score ind]=max(score);
        else
            [score(1),match_only1(1),breakpoint2(1),breakpoint1(1),onread(1),skip(1)]=align_with_breakpoint(seq2, seq1, reads{i}, qual{i}, params.align_wsplit);            
            % reads{i}(1:onread(1)-skip(1))
            % seq2(breakpoint2(1)-onread(1)+1+skip(1):breakpoint2(1))
            % reads{i}(onread(1)+1-skip(1):onread(1))
            % reads{i}(onread(1)+1:end)
            % seq1(breakpoint1(1)+1:breakpoint1(1)+length(reads{i})-onread(1))          
            [score(2),match_only1(2),breakpoint2(2),breakpoint1(2),onread(2),skip(2)]=align_with_breakpoint(seq2, seq1, revreads{i}, qual{i}(end:-1:1), params.align_wsplit);
            [score ind]=max(score);
            % revreads{i}(1:onread(ind)-skip(ind))
            % seq2(breakpoint2(ind)-onread(ind)+1+skip(ind):breakpoint2(ind))
            % revreads{i}(onread(ind)+1:end)
            % seq1(breakpoint1(ind)+1:breakpoint1(ind)+length(reads{i})-onread(ind))
        end
        if ((~match_only1(ind))&&(score>15))
            l=min(sum([p(1:(onread(ind)-skip(ind))) p(onread(ind)+1:end)]),0.66*(length(seq1)+length(seq2)))-params.align_wsplit;
            score=score/l;
            if (score> params.align_min_qual)
                %if ind==1
                %    reads{i}(1:onread(1)-skip(1))
                %    seq1(breakpoint1(1)-onread(1)+1+skip(1):breakpoint1(1))
                %    reads{i}(onread(1)+1-skip(1):end)
                %    seq2(breakpoint2(1)+1:breakpoint2(1)+length(reads{i})-onread(1))
                %end
                %	[breakpoint1(ind) breakpoint2(ind) skip(ind) onread(ind) ind]     
                if (skip(ind)==0)  % when there's no insertion of non templated sequence, breakpoint is not well defined around the microhomoogy (if any), and so always define the rightmost for standartization (after all we do a majorty vote...)
                    len=min(ls2-breakpoint2(ind),ls1-breakpoint1(ind))-1;
                    j2=0;
                    while (j2<len)&&(seq2(breakpoint2(ind)+1+j2)==seq1(breakpoint1(ind)+1+j2))
                        j2=j2+1;
                    end
                    breakpoint1(ind)=breakpoint1(ind)+j2;
                    breakpoint2(ind)=breakpoint2(ind)+j2;
                end
                breakpoint=([breakpoint1(ind);breakpoint2(ind)])+startongenome+seq_start-1;%+[0;skip(ind)];
                if inversion
                    breakpoint(2)=startongenome(2)+seq_end(2)-breakpoint2(ind);%+skip(ind);                    
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
                if (length(scores)>=params.align_enough_reads)
                    [temp1,temp2,n]=unique(breakpoints,'rows');
                    counts=accumarray(n,1);
                    if max(counts)>=params.align_enough_reads
                        break;
                    end
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
        [u2,m2,n2]=unique(forseqs(n==ind(bk)));
        counts2=accumarray(n2,1);
        [t,i2]=max(counts2);
        forseq=u2{i2};
    end
end
