% this is just an oncopanel formatted maf
m=load_struct('/Users/amaro/Downloads/OncopanelSNVs.txt');
% oncopanel specific fixes
samples=unique(m.aggregate_name);
m.tumor_fraction=str2double(m.tumor_fraction);
m.coverage=str2double(m.coverage);
m.alt_count=round(m.tumor_fraction.*m.coverage);
m.ref_count=m.coverage-m.alt_count;
m.position=str2double(m.position);

% probably wanna take this out :)
VHL_limit_up = 60000000;
VHL_limit_down = 1;
%chr3p_lim=60000000;

X.afbeta=zeros(length(samples),101);
% runs for all samples in the maf
for i=1:length(samples)

X.sample{i,1}=samples{i};

% probably want to remove this or iterate over each chrom arm on its own 
current=reorder_struct(m,ismember(m.aggregate_name,samples{i}));
chr3=reorder_struct(current,ismember(current.chromosome,'3')&current.position<VHL_limit_up&current.position>VHL_limit_down);
af=[0:.01:1];

%tracking number of muts used
X.nMut(i,1)=slength(chr3);

for j=1:slength(chr3)
    %get alt and ref for each mutation 
    x2=betapdf(af,chr3.alt_count(j)+1,chr3.ref_count(j)+1); 
    X.afbeta(i,:)=X.afbeta(i,:)+(x2/sum(x2));
    %add them up and look for peaks (in cn consistent reigon (3p) )
end


% finding peaks  
x1=X.afbeta(i,:);
j=2:(length(X.afbeta(i,:))-1);
% the max peak not matching the homozygous mutations is essentially 
kp=af(find((x1(j-1)<=x1(j))&(x1(j+1)<x1(j))));
kp=kp(kp>.20&kp<.80);


if length(kp)>1 && max(kp)>.51
    
   
   kv=af(find((x1(j-1)>=x1(j))&(x1(j+1)>x1(j))));
   p=max(kp);
   vl=max(kv(find(kv<p)));
   vh=(kv(find(kv>p,1,'first')));
% calculating purity by CN LOH or by a deletion
   X.loh_peak(i,1)=p;
  X.pur_high(i,1)=(p-.5)*2;
  X.pur_low(i,1)=1-((1/p)-1);
  X.pur(i,1)=((X.pur_high(i)+X.pur_low(i))/2);
   
else
   X.loh_peak(i,1)=max(kp);
    X.pur(i,1)=0;
    X.pur_low(i,1)=0;
    X.pur_high(i,1)=0;
end

end
figure()
for i=1:slength(X)
   
subplot(6,4,i)
hold on
plot(af,X.afbeta(i,:)./100)
xlim([.2 .8]);
title(X.sample{i});
pk=(X.loh_peak(i));

y_l=ylim;
if pk<.51
plot([pk pk],[0 y_l(2)],'k--');
else
    plot([pk pk],[0 y_l(2)],'r--');
    pk1=.5-(pk-.5);
    plot([pk1 pk1],[0 y_l(2)],'b--');
end
hold off
 end

%i=7 whole chromosome
%i=14



% af=[0:.01:1];
% for i=1:slength(m)
%     %get alt and ref for each mutation %exclude homozygous
%     beta_dist(:,i)=betapdf(af,alt+1,ref+1); %normalize
%     %add them up and look for peaks (in cn consistent reigon (3p) )
% end