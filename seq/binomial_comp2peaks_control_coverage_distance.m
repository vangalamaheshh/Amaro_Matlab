% annotpos,chrinds contain the chromatin marks data. given n different chipseq experiments, annotpos is a cell of size n*1, and chrinds is a cell of size n*23 so that annotpos{i}{(chrinds{i,j}) contains all the peaks loci on chromosome j in experiment i.
% allbkpts is a struct containing fields Tumor_Sample,chr1,pos1,chr2,pos2
% samples contains the samples names (to match allbkpts.Tumor_Sample)
% avcov is the reference coverage (that is also manually mapped in /xchip/cga1/ydrier/assemble/prostate/coverage_average_chr*)
%
%   Yotam Drier 2010

function [scores1,scores2,rate1,bkrate1,rate2,bkrate2,dist,bkdist,actualdist]=binomial_comp2peaks_control_coverage_distance(annotpos,chrinds,allbkpts,samples,avcov,trials)

chrlen=[247249719 242951149 199501827 191273063 180857866 170899992 158821424 146274826 140273252 135374737 134452384 132349534 114142980 106368585 100338915 88827254 78774742 76117153 63811651 62435964 46944323 49691432 154913754 57772954];
dists=[1e3 2.5e3 5e3 7.5e3 1e4 2.5e4 5e4 7.5e4 1e5 2.5e5 5e5 7.5e5 1e6 5e6];
ld=length(dists);
if ~exist('trials','var')
  trials=1000;
end
m=length(samples);
na=length(annotpos);
pdep1=zeros(ld,m,na);
penr1=zeros(ld,m,na);
rate1=zeros(ld,m,na);
bkrate1=zeros(ld,m,na);
pdep2=zeros(ld,m,na);
penr2=zeros(ld,m,na);
rate2=zeros(ld,m,na);
bkrate2=zeros(ld,m,na);
dist=cell(m,na);
bkdist=cell(m,na);
nb=length(allbkpts.chr1);
map=cell(23,11);
bkcov=zeros(nb,2);
actualdist=zeros(nb,na,trials);

fprintf('Processing coverage...');

for i=1:23,fprintf(' %d/%d',i,23);
    load(['/xchip/cga1/ydrier/assemble/prostate/coverage_average_chr' num2str(i) '.mat'],'mapcov')
    for j=1:11
        map{i,j}=find(mapcov==j);
    end
     ci=find(allbkpts.chr1==i);
     for j=1:length(ci);      
         bkcov(ci(j),1)=mapcov(allbkpts.pos1(ci(j)));
     end
     ci=find(allbkpts.chr2==i);
     for j=1:length(ci);          
         bkcov(ci(j),2)=mapcov(allbkpts.pos2(ci(j)));
     end
end, fprintf('\n');

clear mapcov

fprintf('Processing samples...\n');

for s=1:m
    bksubseti=strmatch(samples{s},allbkpts.Tumor_Sample);
    bksubseti=bksubseti(max(allbkpts.chr1(bksubseti),allbkpts.chr2(bksubseti))<24);
    bkpts=structfun(@(x) x(bksubseti),allbkpts,'UniformOutput',false);
    n=length(bksubseti);
    intra=bkpts.chr1==bkpts.chr2;    
    rearr_dist=abs(bkpts.pos2-bkpts.pos1);
    disp(['**** ' samples{s} ' ****'])    
    for a=1:na
        dist{s,a}=dist2annot3(bkpts,annotpos{a},chrinds(a,:));
        bkdist{s,a}=zeros(n,2,trials);
        rand_bkpts.chr1=bkpts.chr1;
        rand_bkpts.chr2=bkpts.chr2;
        for t=1:trials
            rand_bkpts.pos1=zeros(n,1);
            rand_bkpts.pos2=zeros(n,1);
            for i=1:n
                ch=rand_bkpts.chr1(i);
                rand_bkpts.pos1(i)=map{ch,bkcov(bksubseti(i),1)}(ceil(rand()*length(map{ch,bkcov(bksubseti(i),1)})));
                if intra(i)
                    p=mod(rand_bkpts.pos1(i)+[rearr_dist(i);-rearr_dist(i)]-1,chrlen(ch))+1;
                    [cov,w]=max(avcov{ch}.coverage(p));
                    j=0;
                    while (cov<avcov{ch}.coverage(bkpts.pos2(i))-5)&&(j<rearr_dist(i))
                        j=j+100;
                        if mod(j,3)==0
                            j=j+1000;
                        end
                        p=mod(rand_bkpts.pos1(i)+[rearr_dist(i)+j;rearr_dist(i)-j;-rearr_dist(i)+j;-rearr_dist(i)-j]-1,chrlen(ch))+1;
                        [cov,w]=max(avcov{ch}.coverage(p));
                    end
                    while (cov<avcov{ch}.coverage(bkpts.pos2(i))-5)
                        j=j+200;
                        if mod(j,3)==0
                            j=j+3000;
                        end
                        p=mod(rand_bkpts.pos1(i)+[rearr_dist(i)+j;rearr_dist(i)-j;-rearr_dist(i)+j;-rearr_dist(i)-j]-1,chrlen(ch))+1;
                        [cov,w]=max(avcov{ch}.coverage(p));
                    end
                    rand_bkpts.pos2(i)=p(w);
                else
                    rand_bkpts.pos2(i)=map{rand_bkpts.chr2(i),bkcov(bksubseti(i),2)}(ceil(rand()*length(map{rand_bkpts.chr2(i),bkcov(bksubseti(i),2)})));
                end
                actualdist(bksubseti(i),a,t)=abs(rand_bkpts.pos2(i)-rand_bkpts.pos1(i));
            end
            bkdist{s,a}(:,:,t)=dist2annot3(rand_bkpts,annotpos{a},chrinds(a,:));
        end
        for i=1:ld;     
            % one is enough
            p=mean(bkdist{s,a}(:)<dists(i));
            k=sum(dist{s,a}(:)<dists(i));
            pdep1(i,s,a) = binocdf(k,2*n,p);
            penr1(i,s,a) = 1-binocdf(k-1,2*n,p);
            rate1(i,s,a)=k/(2*n);
            bkrate1(i,s,a)=p;
            % require both
            p=mean(mean((bkdist{s,a}(:,1,:)<dists(i))&(bkdist{s,a}(:,2,:)<dists(i)),3)); 
            k=sum(sum((dist{s,a}(:,1,:)<dists(i))&(dist{s,a}(:,2,:)<dists(i)),3));
            pdep2(i,s,a) = binocdf(k,n,p);
            penr2(i,s,a) = 1-binocdf(k-1,n,p);
            rate2(i,s,a)=k/n;
            bkrate2(i,s,a)=p;
        end
    end
end
enr=penr1<pdep1;
p=pdep1;
p(enr)=penr1(enr);
scores1=log2(p);
scores1(enr)=-scores1(enr);
enr=penr2<pdep2;
p=pdep2;
p(enr)=penr2(enr);
scores2=log2(p);
scores2(enr)=-scores2(enr);
