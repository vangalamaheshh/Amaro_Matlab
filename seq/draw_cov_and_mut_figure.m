function draw_cov_and_mut_figure(TCGA)

% make coverage graph

N_cov = squeeze(TCGA.N_cov(:,TCGA.TOT,:));
N_terr = squeeze(TCGA.N_terr(:,TCGA.TOT));
n_tot = squeeze(TCGA.n_nonsilent(:,TCGA.TOT,:));

F = N_cov ./ repmat(N_terr,1,TCGA.np);
gene = TCGA.gene;
patient = TCGA.patient;
if ~isnumeric(gene.phase), gene.phase = str2double(gene.phase); end

% sort by center
[tmp ord] = sort(gene.center);
F = F(ord,:);
n_tot = n_tot(ord,:);
gene = reorder_struct(gene,ord);

% sort by phase
[tmp ord] = sort(gene.phase);
F = F(ord,:);
n_tot = n_tot(ord,:);
gene = reorder_struct(gene,ord);

% sort by whether patients were treated, then by whether in first batch of 100
tr = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/treated_samples.txt');
tmp = ismember(patient.name,tr.name);
ord = [find(tmp);find(~tmp)];
F = F(:,ord);
n_tot = n_tot(:,ord);
patient = reorder_struct(patient,ord);
p100 = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/100patients.txt');
tmp = ismember(patient.name,p100.name);
ord = [find(tmp);find(~tmp)];
F = F(:,ord);
n_tot = n_tot(:,ord);
patient = reorder_struct(patient,ord);

% display coverage graph

close; clf; hold on; imagesc(F); colorbar;
y = []; n = {};
for g=1:TCGA.ng
  if g==1 || ~strcmp(gene.center{g},gene.center{g-1}) ||...
             gene.phase(g)~=gene.phase(g-1)
    y = [y g];
    label = sprintf('Phase%d/%s',gene.phase(g),gene.center{g});
    n = [n; label];
  end
end
set(gca,'YTick',y);
set(gca,'YTickLabel',n);
set(gca, 'ydir','reverse')
set(gca, 'ylim',[0 TCGA.ng]);
set(gca, 'xlim',[0 TCGA.np]);
xlabel('patient');
ylabel('gene');
set(gca,'position',[0.15 0.1 0.75 0.85]);
set(gcf,'position',[103 21 1136 915]);
set(gcf,'color',[1 1 1]);
% superimpose mutations
for p=1:TCGA.np, for g=1:TCGA.ng
  if n_tot(g,p)
    rectangle('position',[p-0.5 g-0.5 1 10],'curvature',[1 1],'facecolor',[0 0 0]);
    rectangle('position',[p-0.5 g-0.5 0.8 8],'curvature',[1 1],'facecolor',[1 1 1]);
  end
end, end
hold off



