call=load_struct('/Users/amaro/Downloads/SH4885Tumor.cov');
m.alt_count=call.i_t_alt_count;
m.alt_count=str2double(call.i_t_alt_count);
m.ref_count=str2double(call.i_t_ref_count);
m.position=str2double(call.Start_position);

m.chromosome=call.Chromosome;
current=m;
i=9;
af=[0:.01:1];

VHL_limit_up = 1200000000;
VHL_limit_down = 1;
%%VHL_limit_up = 20197354;
%VHL_limit_down = 9181319;
chr3=reorder_struct(current,ismember(current.chromosome,'3')&current.position<VHL_limit_up&current.position>VHL_limit_down);
X.nMut(i,1)=slength(chr3);


X.Sample{i,1}='4885';

X.nMut(i,1)=slength(chr3);
for j=1:slength(chr3)
%get alt and ref for each mutation %exclude homozygous
x2=betapdf(af,chr3.alt_count(j)+1,chr3.ref_count(j)+1); %normalize
X.afbeta(i,:)=X.afbeta(i,:)+(x2/sum(x2));
%add them up and look for peaks (in cn consistent reigon (3p) )
end
x1=X.afbeta(i,:);
j=2:(length(X.afbeta(i,:))-1);
kp=af(find((x1(j-1)<=x1(j))&(x1(j+1)<x1(j))));
kp=kp(kp>.20&kp<.80);


kv=af(find((x1(j-1)>=x1(j))&(x1(j+1)>x1(j))));
p=max(kp);
vl=max(kv(find(kv<p)));
vh=(kv(find(kv>p,1,'first')));
%   p_dist=cumsum(x1(find(af==vl):find(af==vh)))./sum(x1(find(af==vl):find(af==vh)));
%   af_h=find(p_dist<.975,1,'last');
%   af_l=find(p_dist>.025,1,'first');
%   p_l=vl+af(af_l);
%   p_h=vl+af(af_h);
X.loh_peak(i,1)=p;
X.pur_high(i,1)=(p-.5)*2;
X.pur_low(i,1)=1-((1/p)-1);
X.pur(i,1)=((X.pur_high(i)+X.pur_low(i))/2);