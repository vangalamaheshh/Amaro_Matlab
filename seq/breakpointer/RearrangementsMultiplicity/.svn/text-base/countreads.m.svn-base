function support=countreads(dRangerFile,BreakPointerFile,bamfile,startfrom,endat,P) 
% Yotam Drier, yotamd@gmail.com

%drangerFile='~/xchip/assemble/countreads/mm/MM-0392/MM-0392-Tumor.breakpoints.txt';
%clear java
%javaclasspath('/xchip/cga1/ydrier/assemble/GrabSplitReads/BamGrasp/GrabSplitReads.jar')
%bamfile='/xchip/cga1/firehose_output/breakpointer_paper/Individual/MM-0392/wgs/bam/raw/tumor.bam';

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'quiet',1);
P=impose_default_value(P,'maxreads',1e7);
P=impose_default_value(P,'blacklist','none');
P=impose_default_value(P,'pairs_extraction',1000); 
P=impose_default_value(P,'readlen',101); 
P=impose_default_value(P,'gsr_jar','/xchip/cga1/ydrier/CancerGenomeAnalysis/trunk/analysis_pipeline/tools/dist/GrabSplitReads.jar');
%P=impose_default_value(P,'gsr_jar','/xchip/cga1/ydrier/assemble/GrabSplitReads/BamGrasp/GrabSplitReads.jar');
P=impose_default_value(P,'joindump_cmd','/xchip/cga1/ydrier/assemble/sandbox/ydrier/breakpointer_analysis/code/RearrangementsMultiplicity/joindump.pl');
P=impose_default_value(P,'refdir','/xchip/tcga/gbm/analysis/lawrence/genome/hg18');
P=impose_default_value(P,'align_min_qual',0.75);
P=impose_default_value(P,'minscorediff',20);
P=impose_default_value(P,'min2span',8);

javaclasspath(P.gsr_jar);
import org.broadinstitute.cga.tools.*;
import java.lang.*;
try gsr = org.broadinstitute.cga.tools.seq.GrabSplitReads(String(bamfile),P.blacklist,P.maxreads,String(P.refdir));
catch gsr = org.broadinstitute.cga.tools.seq.GrabSplitReads(String(bamfile),P.blacklist,P.maxreads,String(P.refdir)); %#ok<CTCH>
end

flds1={'num','somatic_score'};
%flds1{1}='drnum'; % dRangerFile can also be BreakPointer output, if filtering by somatic_score is not needed
dr=make_numeric(keep_fields(load_struct(dRangerFile),flds1),flds1);
%dr.num=dr.drnum; dr = rmfield(dr, 'drnum');
flds2={'drnum','chr1','chr2','pos1','pos2','kept1','kept2','fs'};
bp=make_numeric(keep_fields(load_struct(BreakPointerFile),flds2),{'drnum','chr1','chr2','pos1','pos2'});
bpsuccess=find(~strcmp(bp.kept1,'failed'));
drsomatic=find(dr.somatic_score>0);
[c,id,ib]=intersect(dr.num(drsomatic),bp.drnum(bpsuccess));
rearr=merge_structs2({structfun(@(x)x(drsomatic(id)),dr,'UniformOutput',false);structfun(@(x)x(bpsuccess(ib)),bp,'UniformOutput',false)});
if ~exist('startfrom','var'), startfrom = 1; end
if ~exist('endat','var'), endat = length(c); end
if endat==-1 || endat>length(c), endat = length(c); end
if startfrom>endat
    support=[];
else
    warning('OFF','Bioinfo:swalign:EmptyAlignment');
    l=endat-startfrom+1;
    support = struct('drnum', rearr.drnum(startfrom:endat), 'first_by_pair', ...
                     zeros(l,1), 'second_by_pair', zeros(l,1), 'fusion_by_pair', ...
                     zeros(l,1), 'first_by_span', zeros(l,1), 'second_by_span', ...
                     zeros(l,1), 'fusion_by_span', zeros(l,1), 'contradict', ...
                     zeros(l,1), 'undecided', zeros(l,1), 'error', zeros(l,1));    
    for i=startfrom:endat
        dname = num2str(i);
        if ~exist(dname,'dir'), mkdir(dname); end
        start1=rearr.pos1(i)-P.pairs_extraction;
        end1=rearr.pos1(i)+P.pairs_extraction;       
        if ~exist([dname '/first.txt'],'file')
            gsr.WritePairs(rearr.chr1(i),start1,end1,[dname '/first.txt']);
        end
        start2=rearr.pos2(i)-P.pairs_extraction;
        end2=rearr.pos2(i)+P.pairs_extraction;        
        if ~exist([dname '/second.txt'],'file')          
            gsr.WritePairs(rearr.chr2(i),start2,end2,[dname '/second.txt']);
        end      
        if ~exist([dname '/bothends.txt'],'file')
            unix(['cat ' dname '/first.txt ' dname '/second.txt | sort | perl ' P.joindump_cmd ' 1 > ' dname '/bothends.txt']);
        end
        pairs=load_struct_noheader([dname '/bothends.txt'],15,{'readname';'chr1';'start1';'end1';'strand1';'qual1';'seq1';'bqual1';'chr2';'start2';'end2';'strand2';'qual2';'seq2';'bqual2'});
        pairs=make_numeric(pairs,{'chr1';'start1';'end1';'strand1';'chr2';'start2';'end2';'strand2'});
        fused=[rearr.kept1{i} rearr.fs{i} rearr.kept2{i}];
        ref1all=genome_region_fun(rearr.chr1(i),start1-P.readlen,end1+P.readlen,P.refdir);
        ref2all=genome_region_fun(rearr.chr2(i),start2-P.readlen,end2+P.readlen,P.refdir);
        rev1all=seqrcomplement(ref1all);
        rev2all=seqrcomplement(ref2all);
        revf=seqrcomplement(fused);
        refs={ref1all;rev1all;ref2all;rev2all;fused;revf};
        l=length(pairs.seq1);
        k=i-startfrom+1;
        for j=1:l
            [first,half1]=classpair(pairs.seq1{j},refs,P.minscorediff,P.align_min_qual,P.min2span);
            [second,half2]=classpair(pairs.seq2{j},refs,P.minscorediff,P.align_min_qual,P.min2span);
            bothsides=(((half1==0)&&(half2==1))||((half1==1)&&(half2==0)));
            if (first==4)&&(second==4)&&bothsides
                support.first_by_pair(k)=support.first_by_pair(k)+1;
                if (~P.quiet), fprintf('pair1: %s (j=%d)\n',pairs.readname{j},j); end
            else
                if (first==5)&&(second==5)&&bothsides
                    support.second_by_pair(k)=support.second_by_pair(k)+1;
                    if (~P.quiet), fprintf('pair2: %s (j=%d)\n',pairs.readname{j},j); end
                else
                    if ((first==4)&&(second==5))||((first==5)&&(second==4))
                        if (half1==2)||(half2==2)
                            support.contradict(k)=support.contradict(k)+1;
                            if (~P.quiet), fprintf('contradict 1: %s (j=%d)\n', pairs.readname{j},j); end
                        else
                            support.fusion_by_pair(k)=support.fusion_by_pair(k)+1;
                            if (~P.quiet), fprintf('pair3: %s (j=%d)\n', pairs.readname{j},j); end
                        end
                    end
                end
            end
            if ((first==3)||(second==3))
                support.fusion_by_span(k)=support.fusion_by_span(k)+1;
                if (~P.quiet), fprintf('single3: %s (j=%d)\n', pairs.readname{j},j); end
            else
                if ((first==6)||(second==6))
                    support.error(k)=support.error(k)+1;
                    if (~P.quiet), fprintf('error: %s (j=%d)\n', pairs.readname{j},j); end
                else
                    if (first==1)||(second==1)
                        if (first==2)||(first==5)||(second==2)||(second==5)
                            support.contradict(k)=support.contradict(k)+1;
                            if (~P.quiet), fprintf('contradict 2: %s (j=%d)\n', pairs.readname{j},j); end
                        else
                            support.first_by_span(k)= support.first_by_span(k)+1;
                            if (~P.quiet), fprintf('single1: %s (j=%d)\n', pairs.readname{j},j); end
                        end
                    else
                        if (first==2)||(second==2)
                            if (first==1)||(first==4)||(second==1)||(second==4)
                                support.contradict(k)=support.contradict(k)+1;
                                if (~P.quiet), fprintf('contradict 3: %s (j=%d)\n', pairs.readname{j},j); end
                            else
                                support.second_by_span(k)=support.second_by_span(k)+1;
                                if (~P.quiet), fprintf('single2: %s (j=%d)\n',pairs.readname{j},j); end
                            end
                        else
                            if (first<=-2)||(second<=-2)
                                support.undecided(k)=support.undecided(k)+1;
                                if (~P.quiet), fprintf('undecided: %s (j=%d)\n',pairs.readname{j},j); end
                            end
                        end
                    end
                end
            end
        end
    end
    warning('ON','Bioinfo:swalign:EmptyAlignment');
end
