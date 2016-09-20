function [x X] = calculate_x_and_X_tables(M,G)



% data about category-specific mutation rate
x_c = sum(sum(M.nnon+M.nsil+M.nflank,3),1);
X_c = sum(sum(M.Nnon+M.Nsil+M.Nflank,3),1);
mu_c = (x_c./X_c);
x_tot = sum(x_c);
X_tot = X_c(end);
Xrel_c = X_c/X_tot;
mu_tot = x_tot/X_tot;
murel_c = mu_c/mu_tot;

% data about sample-specific mutation rate
x_s = sum(sum(M.nnon+M.nsil+M.nflank,2),1);
X_s = sum(M.Nnon(:,end,:)+M.Nsil(:,end,:)+M.Nflank(:,end,:),1);
Xrel_s = X_s/sum(X_s);
mu_s = (x_s./X_s);
mu_tot = (sum(x_s)/sum(X_s));
murel_s = mu_s/mu_tot;

% as a CHECK, see how well the marginals predict the joint distribution
x_cs = sum(M.nnon,1);
X_cs = sum(M.Nnon,1);
mu_cs_obs = x_cs ./ X_cs;
murel_cs_obs = mu_cs_obs / mu_tot;
murel_cs_fit = bsxfun(@times,murel_c,murel_s);
catno = repmat(1:8,[1 1 size(M.nnon,3)]);
x = log10(murel_cs_obs(:)); x(murel_cs_obs(:)==0) = -3.5;
y = log10(murel_cs_fit(:)); y(murel_cs_fit(:)==0) = -3.5;
clr = catno(:);
clf;scatter(x,y,[],clr);line(xlim,xlim); % ---> looks good!


% data about gene-specific mutation rate
x_gene = G.nsil + G.nflank;
X_gene = G.Nsil + G.Nflank;
x_bagel = G.nfit - x_gene;
X_bagel = G.Nfit - X_gene;

% combine these into a single x+X:
w_bagel = ones(length(x_gene),1); w_bagel(G.Fmle>1.4) = 0.2;
x_sphere = x_gene + w_bagel.*x_bagel;
X_sphere = X_gene + w_bagel.*X_bagel;

% add prior from what the per-gene rates are like
[amle,bmle]=mle_beta(x_gene,X_gene);    % 1.72, 251950
x_sphere = x_sphere + (amle-1);
X_sphere = X_sphere + (bmle+amle-2);

%% scale by mu_cs
x_sphere_cs = bsxfun(@times,x_sphere,murel_cs_fit);
X_sphere_cs = repmat(X_sphere,[1 size(x_sphere_cs,2) size(x_sphere_cs,3)]);

x = x_sphere_cs;
X = X_sphere_cs;

