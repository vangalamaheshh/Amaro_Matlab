function [ge,cn,ps,snps]=gene_summary(Cr,Gr,rkg,ext,gsymb,greg,bad_probes,use_flanking)

ge=[];
cn=[];
ps=[];
snps=[];

if ~exist('use_flanking','var')
  use_flanking=0;
end

if ~ischar(rkg)
  ki=grep(['^' gsymb],rkg(:,1),1);
else
  ki=[];
end

% [Cg,v]=collapse_to_genes(Cr,rg,rkg{ki,4});
ps=grep(['^' gsymb '$'],Gr.gsymb,1);

if ~exist('bad_probes','var') 
  figure(1);clf
  subplot(2,1,1)
  imagesc(Gr.dat(ps,:));
  subplot(2,1,2)
  imagesc(ctr_norm(Gr.dat(ps,:),'median','mad'));
  return
elseif ~isempty(bad_probes)
  ps(bad_probes)=[]; % bad probes
end

if ~isstruct(greg)
  peak=greg;
else
  peak=greg.peak;
end

if ~isempty(ki)
  snps=find_snps(Cr,rkg{ki,6},[],[],use_flanking);
else
  snps=find_snps(Cr,rkg,[],[],use_flanking);
end

if peak<0
  peak=floor(mean(snps));
end

cn=Cr.dat(peak,:);
ge=median(ctr_norm(Gr.dat(ps,:),'median','mad'),1);
fs=14;  
figure(1); clf;
plot(cn,ge,'b.'); hold on
plot(mean(Cr.dat(snps,:),1),ge,'ro');
title(gsymb,'FontSize',fs);
% legend({'gene location'},'Location','northwest');
print_D([gsymb '_cn_ge.mm' ext],{{'pdf'},{'fig'}});










