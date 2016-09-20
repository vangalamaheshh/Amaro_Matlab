function S=gcs2n(dat,cls0,cls1)

m0=nanmean(dat(:,cls0),2);
s0=nanstd(dat(:,cls0)')';
m1=nanmean(dat(:,cls1),2);
s1=nanstd(dat(:,cls1)')';
S=(m1-m0)./(gc_std(m0,s0)+gc_std(m1,s1));

