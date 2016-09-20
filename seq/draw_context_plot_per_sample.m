function draw_context_plot(patname,R,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'fontsize',20);

if length(patname)~=length(R), error('length(patname)~=length(R)'); end

[tmp ui] = unique(R{1}.base_context); bc = R{1}.base_context(sort(ui));
rate = nan(length(patname),length(bc),3);
for i=1:length(patname)
  r = reorder_struct(R{i},strcmp('total',R{i}.regulatory_potential) & strcmp('total',R{i}.transcript_zone));
  for j=1:length(bc)
    idx = find(strcmp(bc{j},r.base_context));
    rate(i,j,1) = r.mutrate(idx);
    rate(i,j,2) = r.ci_low(idx);
    rate(i,j,3) = r.ci_high(idx);
end,end
% sort by total rate
[tmp ord] = sort(rate(:,end,1));
rate = rate(ord,:,:);
patname = patname(ord);
bc = regexprep(bc,'_',' ');
% draw figure
figure(1),clf, fs = P.fontsize;
subplot('position',[0.1 0.1 0.12 0.8]);
barweb(rate(:,1,1)',rate(:,1,3)',rate(:,1,2)',0.8,bc(1)); ylabel('mutations per million sites','fontsize',fs);
set(gca,'fontsize',fs);
subplot('position',[0.26 0.1 0.67 0.8]);
barweb(rate(:,2:end,1)',rate(:,2:end,3)',rate(:,2:end,2)',0.8,bc(2:end)); ylabel('mutations per million sites', ...
                                                  'fontsize',fs);
set(gca,'fontsize',fs);
set(gcf,'position',[1000 20 2399 904],'color',[1 1 1]);
legend(patname,'location','best');

