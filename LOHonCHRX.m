% 
% 
% for j=1:slength(chr3)
% %get alt and ref for each mutation %exclude homozygous
% x2=betapdf(af,chr3.alt_count(j)+1,chr3.ref_count(j)+1); %normalize
% X.afbeta(i,:)=X.afbeta(i,:)+(x2/sum(x2));
% %add them up and look for peaks (in cn consistent reigon (3p) )
% end
% x1=X.afbeta(i,:);
% j=2:(length(X.afbeta(i,:))-1);
% kp=af(find((x1(j-1)<=x1(j))&(x1(j+1)<x1(j))));
% kp=kp(kp>.20&kp<.80);
% 
% 
% kv=af(find((x1(j-1)>=x1(j))&(x1(j+1)>x1(j))));
% p=max(kp);
% vl=max(kv(find(kv<p)));
% vh=(kv(find(kv>p,1,'first')));
% %   p_dist=cumsum(x1(find(af==vl):find(af==vh)))./sum(x1(find(af==vl):find(af==vh)));
% %   af_h=find(p_dist<.975,1,'last');
% %   af_l=find(p_dist>.025,1,'first');
% %   p_l=vl+af(af_l);
% %   p_h=vl+af(af_h);
% X.loh_peak(i,1)=p;
% X.pur_high(i,1)=(p-.5)*2;
% X.pur_low(i,1)=1-((1/p)-1);
% X.pur(i,1)=((X.pur_high(i)+X.pur_low(i))/2);


gender=load_struct('/Users/amaro/Documents/CLL_Questions_For_Dan/XchrData_gender_assignment.txt');
af=[0:.01:1];

for i=1:slength(gender)
    if strcmp(gender.gender{i},'Female')
        m=load_struct(gender.germline_maf_singlesample_analysisready{i});
        m=reorder_struct(m,~ismember(m.allelic_depth,''));
        m=reorder_struct(m,ismember(m.Chromosome,'X'));       
        for j=1:slength(m)
        strs=split(m.allelic_depth{j},',');
        m.alt_count(j,1)=str2double(strs{2});
        m.ref_count(j,1)=str2double(strs{1});
        x2=betapdf(af,m.alt_count(j)+1,m.ref_count(j)+1);
        m.afbeta(j,:)=(x2/sum(x2));
        end
        m.af=(m.alt_count./(m.alt_count+m.ref_count));
        m.pos=xhg19(m.Chromosome,str2double(m.Start_position));
        j=2:(length(m.afbeta(i,:))-1);
        x1=(sum(m.afbeta));
        kp=af(find((x1(j-1)<=x1(j))&(x1(j+1)<x1(j))));
        gender.LOHstat(i,1)=min(abs(kp-.5));
        i
    end
end