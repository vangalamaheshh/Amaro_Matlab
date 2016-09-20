function V = load_vog()
% load all Vogelstein data

A = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');

tmpdir = pwd;
cd /xchip/tcga/gbm/analysis/lawrence/vog;

% breast and colorectal

fprintf('Loading breast and colorectal data\n');
V.brco.targ = load_struct('brco/tableS1.txt');
% remove intrusive hyphens from some of the target gene names
V.brco.targ.gene = regexprep(V.brco.targ.gene,'-$','');
V.brco.targ.gene = regexprep(V.brco.targ.gene,'^TMEM13-2B$','TMEM132B');
V.brco.targ.gene = regexprep(V.brco.targ.gene,'^FAM13-A1$','FAM13A1');
V.brco.targ.gene = regexprep(V.brco.targ.gene,'^FAM13-C1$','FAM13C1');
V.brco.targ.gene = apply_aliases(V.brco.targ.gene,A);
V.brco.targ.chr = convert_chr(regexprep(V.brco.targ.site,':.*',''));
tmp = regexprep(V.brco.targ.site,'.*:','');
V.brco.targ.start = str2double(regexprep(tmp,'-.*',''));
V.brco.targ.end = str2double(regexprep(tmp,'.*-',''));
%%%%%%%% NOTE: targets are still on hg17 !

V.brco.mut_discval = load_struct('brco/tableS3.txt');
V.brco.mut_prev = load_struct('brco/tableS5.txt','%s%s%f%f%f%f%f');
V.brco.mut_discval.gene = apply_aliases(V.brco.mut_discval.gene,A);
V.brco.mut_prev.gene = apply_aliases(V.brco.mut_prev.gene,A);

V.brco.ns_disc = 11;
V.brco.ns_val = 24;
V.brco.ns_prev = 96;

[V.brco.gene.name gi V.brco.targ.gene] = unique(V.brco.targ.gene);
V.brco.ng = slength(V.brco.gene);
V.brco.gene.terr = compute_terr(V.brco);

   function terr = compute_terr(Vx)
     terr = zeros(Vx.ng,1);
     for i=1:Vx.ng
       if ~mod(i,1000), fprintf('%d/%d ',i,Vx.ng); end
       idx = find(Vx.targ.gene==i);
       s = Vx.targ.start(idx);
       e = Vx.targ.end(idx);
       mn = min(s);
       mx = max(e);
       len = mx-mn+1;
       a = zeros(len,1);
       for k=1:length(idx);
         j=idx(k);
         s2 = Vx.targ.start(j);
         e2 = Vx.targ.end(j);
         a(s2-mn+1:e2-mn+1)=1;
       end
       terr(i) = sum(a);
     end
     fprintf('\n');
   end

V.brco.mut_discval.gene = listmap(V.brco.mut_discval.gene, V.brco.gene.name);
V.brco.gene.n_br_disc = zeros(V.brco.ng,1);
V.brco.gene.n_br_val = zeros(V.brco.ng,1);
V.brco.gene.n_co_disc = zeros(V.brco.ng,1);
V.brco.gene.n_co_val = zeros(V.brco.ng,1);
for i=1:slength(V.brco.mut_discval)
  if strcmp(V.brco.mut_discval.type{i},'Synonymous'), continue; end
  g = V.brco.mut_discval.gene(i);
  switch(V.brco.mut_discval.tumor_type{i})
    case 'Breast'
      switch(V.brco.mut_discval.screen{i})
        case 'Discovery'
          V.brco.gene.n_br_disc(g) = V.brco.gene.n_br_disc(g) + 1;
        case 'Validation'
          V.brco.gene.n_br_val(g) = V.brco.gene.n_br_val(g) + 1;
        otherwise
          error('Unknown screen type in V.brco.mut_discval')
      end
    case 'Colorectal'
      switch(V.brco.mut_discval.screen{i})
        case 'Discovery'
          V.brco.gene.n_co_disc(g) = V.brco.gene.n_co_disc(g) + 1;
        case 'Validation'
          V.brco.gene.n_co_val(g) = V.brco.gene.n_co_val(g) + 1;
        otherwise
          error('Unknown screen type in V.brco.mut_discval')
      end
    otherwise
      error('Unknown tumor type in V.brco.mut_discval')
  end
end

V.brco.gene.n_co_val(V.brco.gene.n_co_disc==0) = nan;
V.brco.gene.n_br_val(V.brco.gene.n_br_disc==0) = nan;

V.brco.mut_prev.gene = listmap(V.brco.mut_prev.gene, V.brco.gene.name);
V.brco.gene.n_co_prev = nan(V.brco.ng,1);
V.brco.gene.n_co_prev(V.brco.mut_prev.gene) = V.brco.mut_prev.n_nonsilent_in_96;

V.brco.gene.r_br_disc = V.brco.gene.n_br_disc ./ (11*V.brco.gene.terr);
V.brco.gene.r_br_val = V.brco.gene.n_br_val ./ (24*V.brco.gene.terr);
V.brco.gene.r_co_disc = V.brco.gene.n_co_disc ./ (11*V.brco.gene.terr);
V.brco.gene.r_co_val = V.brco.gene.n_co_val ./ (24*V.brco.gene.terr);
V.brco.gene.r_co_prev = V.brco.gene.n_co_prev ./ (96*V.brco.gene.terr);


% GBM

fprintf('Loading GBM data\n');
V.gbm.targ=load_struct('gbm/tableS1.txt');
V.gbm.targ.chr = convert_chr(regexprep(V.gbm.targ.site,':.*',''));
tmp = regexprep(V.gbm.targ.site,'.*:','');
V.gbm.targ.start = str2double(regexprep(tmp,'-.*',''));
V.gbm.targ.end = str2double(regexprep(tmp,'.*-',''));
V.gbm.targ.gene = apply_aliases(V.gbm.targ.gene,A);

V.gbm.mut_disc = load_struct('gbm/tableS3.txt');
% remove hypermutated tumor
V.gbm.mut_disc = reorder_struct(V.gbm.mut_disc,~strcmpi(V.gbm.mut_disc.tumor,'Br27P'));
V.gbm.mut_disc.gene = apply_aliases(V.gbm.mut_disc.gene,A);

V.gbm.mut_prev = load_struct('gbm/tableS4.txt');
V.gbm.mut_prev.gene = apply_aliases(V.gbm.mut_prev.gene,A);

V.gbm.ns_disc = 22;
V.gbm.ns_prev = 83;

[V.gbm.gene.name gi V.gbm.targ.gene] = unique(V.gbm.targ.gene);
V.gbm.ng = slength(V.gbm.gene);
V.gbm.gene.terr = compute_terr(V.gbm);

V.gbm.mut_disc.gene = listmap(V.gbm.mut_disc.gene,V.gbm.gene.name);
V.gbm.mut_prev.gene = listmap(V.gbm.mut_prev.gene,V.gbm.gene.name);

V.gbm.gene.n_disc = zeros(V.gbm.ng,1);
for i=1:slength(V.gbm.mut_disc)
  if strcmp(V.gbm.mut_disc.type{i},'Synonymous'), continue; end
  g = V.gbm.mut_disc.gene(i);
  V.gbm.gene.n_disc(g) = V.gbm.gene.n_disc(g) + 1;
end

V.gbm.gene.n_prev = nan(V.gbm.ng,1);
for i=1:slength(V.gbm.mut_prev)
  if strcmp(V.gbm.mut_prev.type{i},'Synonymous'), continue; end
  g = V.gbm.mut_prev.gene(i);
  if isnan(V.gbm.gene.n_prev(g)), V.gbm.gene.n_prev(g) = 0; end
  V.gbm.gene.n_prev(g) = V.gbm.gene.n_prev(g) + 1;
end

idx = listmap({'ASTN', 'KIAA0133', 'KIAA0774', 'SERPINA12', 'TRPV5'},V.gbm.gene.name);
V.gbm.gene.n_prev(idx) = 0;   % these five genes were in the prevalence screen but had 0 mutations

V.gbm.gene.r_disc = V.gbm.gene.n_disc ./ (21*V.gbm.gene.terr);
V.gbm.gene.r_prev = V.gbm.gene.n_prev ./ (83*V.gbm.gene.terr);

% pancreatic

fprintf('Loading pancreatic data\n');
V.panc.targ = load_struct('panc/tableS2.txt');
V.panc.targ.chr = convert_chr(regexprep(V.panc.targ.site,':.*',''));
tmp = regexprep(V.panc.targ.site,'.*:','');
V.panc.targ.start = str2double(regexprep(tmp,'-.*',''));
V.panc.targ.end = str2double(regexprep(tmp,'.*-',''));
V.panc.targ.gene = apply_aliases(V.panc.targ.gene,A);

V.panc.mut_disc = load_struct('panc/tableS3.txt');
V.panc.mut_prev = load_struct('panc/tableS4.txt');
V.panc.mut_disc.gene = apply_aliases(V.panc.mut_disc.gene,A);
V.panc.mut_prev.gene = apply_aliases(V.panc.mut_prev.gene,A);

V.panc.ns_disc = 24;
V.panc.ns_prev = 90;

[V.panc.gene.name gi V.panc.targ.gene] = unique(V.panc.targ.gene);
V.panc.ng = slength(V.panc.gene);
V.panc.gene.terr = compute_terr(V.panc);

V.panc.mut_disc.gene = listmap(V.panc.mut_disc.gene,V.panc.gene.name);
V.panc.mut_prev.gene = listmap(V.panc.mut_prev.gene,V.panc.gene.name);

V.panc.gene.n_disc = zeros(V.panc.ng,1);
for i=1:slength(V.panc.mut_disc)
  if strcmp(V.panc.mut_disc.type{i},'Synonymous'), continue; end
  g = V.panc.mut_disc.gene(i);
  V.panc.gene.n_disc(g) = V.panc.gene.n_disc(g) + 1;
end

V.panc.gene.n_prev = nan(V.panc.ng,1);
for i=1:slength(V.panc.mut_prev)
  g = V.panc.mut_prev.gene(i);
  if isnan(V.panc.gene.n_prev(g)), V.panc.gene.n_prev(g) = 0; end
  if strcmp(V.panc.mut_prev.type{i},'Synonymous'), continue; end
  V.panc.gene.n_prev(g) = V.panc.gene.n_prev(g) + 1;
end

idx = listmap({'Q9H5F0_HUMAN', 'DPP6', 'KLHDC4', 'DEPDC2', 'FLJ46481', 'MIZF',...
               'RASSF6', 'TM7SF4', 'XR_017918.1', 'MEP1A', 'CNTN5', 'KIAA1024',...
               'SLC1A6', 'OVCH1', 'ADAMTS20', 'OR10R2'},V.panc.gene.name);
V.panc.gene.n_prev(idx) = 0;   % sequenced in prev screen but had 0 mutations

V.panc.gene.r_disc = V.panc.gene.n_disc ./ (24*V.panc.gene.terr);
V.panc.gene.r_prev = V.panc.gene.n_prev ./ (90*V.panc.gene.terr);

% combine all data

fprintf('Combining all data\n');
V.gene=[];
V.gene.name = union(V.brco.gene.name, [V.gbm.gene.name; V.panc.gene.name]);
V.ng = slength(V.gene);

V.screen.name = {'br_disc';'br_val';'co_disc';'co_val';'co_prev';...
  'gbm_disc';'gbm_prev';'panc_disc';'panc_prev'};
V.nscr = slength(V.screen);
V.screen.ns = [V.brco.ns_disc;V.brco.ns_val;V.brco.ns_disc;V.brco.ns_val;V.brco.ns_prev;...
               V.gbm.ns_disc;V.gbm.ns_prev;V.panc.ns_disc;V.panc.ns_prev];

V.n_nonsilent = nan(V.ng,V.nscr);
V.N_terr = nan(V.ng,V.nscr);

idx = listmap(V.gene.name,V.brco.gene.name);
V.gene.brco_idx = idx;
i1 = find(~isnan(idx));
i2 = idx(i1);
V.n_nonsilent(i1,1) = V.brco.gene.n_br_disc(i2);
V.n_nonsilent(i1,2) = V.brco.gene.n_br_val(i2);
V.n_nonsilent(i1,3) = V.brco.gene.n_co_disc(i2);
V.n_nonsilent(i1,4) = V.brco.gene.n_co_val(i2);
V.n_nonsilent(i1,5) = V.brco.gene.n_co_prev(i2);
V.N_terr(i1,1:5) = repmat(V.brco.gene.terr(i2),1,5);

idx = listmap(V.gene.name,V.gbm.gene.name);
V.gene.gbm_idx = idx;
i1 = find(~isnan(idx));
i2 = idx(i1);
V.n_nonsilent(i1,6) = V.gbm.gene.n_disc(i2);
V.n_nonsilent(i1,7) = V.gbm.gene.n_prev(i2);
V.N_terr(i1,6:7) = repmat(V.gbm.gene.terr(i2),1,2);

idx = listmap(V.gene.name,V.panc.gene.name);
V.gene.panc_idx = idx;
i1 = find(~isnan(idx));
i2 = idx(i1);
V.n_nonsilent(i1,8) = V.panc.gene.n_disc(i2);
V.n_nonsilent(i1,9) = V.panc.gene.n_prev(i2);
V.N_terr(i1,8:9) = repmat(V.panc.gene.terr(i2),1,2);

% multiply territory by number-of-samples (ns) to get coverage

V.N_cov = V.N_terr .* repmat(V.screen.ns',V.ng,1);
V.N_cov(isnan(V.n_nonsilent)) = nan;

% BMRs

V.screen.bmr_excl = nan(9,1);

idx=listmap({'TP53';'PTEN';'RB1'},V.gbm.gene.name);
i2=setdiff(1:V.gbm.ng,idx);
n=sum(V.gbm.gene.n_disc(i2));
N=sum(V.gbm.gene.terr(i2));
r=0.95*n/(21*N);
V.screen.bmr_excl(6:7) = r;

idx=listmap({'SMAD4';'CDKN2A';'TP53';'KRAS'},V.panc.gene.name);
i2=setdiff(1:V.panc.ng,idx);
n=sum(V.panc.gene.n_disc(i2));
N=sum(V.panc.gene.terr(i2));
r=0.96*n/(24*N);
V.screen.bmr_excl(8:9) = r;

V.screen.bmr_lower = [1.40;0.74;0.99;1.44;1.44;0.38;0.38;0.54;0.54]*1e-6;
V.screen.bmr_upper = [3.62;1.91;2.35;3.41;3.41;1.02;1.02;1.38;1.38]*1e-6;
V.screen.bmr_mid = (V.screen.bmr_lower + V.screen.bmr_upper) / 2;

% correlate with TCGA data

fprintf('Correlating with TCGA data\n');
TCGA = load_TCGA;
a = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');
vname = apply_aliases(V.gene.name,a);
idx = listmap(vname,TCGA.gene.name);
i1 = find(~isnan(idx));
i2 = idx(i1);
V.gene.tcga_phase = zeros(V.ng,1);
V.gene.tcga_phase(i1) = TCGA.gene.phase(i2);
V.gene.tcga_phase(isnan(V.gene.tcga_phase))=0;
V.gene.tcga_n_muts = nan(V.ng,1);
V.gene.tcga_n_muts(i1) = sum(TCGA.n_nonsilent(i2,TCGA.TOT,:),3);
V.gene.tcga_n_muts(V.gene.tcga_phase==0) = nan;

% correlate with TSP data

fprintf('Correlating with TSP data\n');
TSP = load_TSP;
a = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');
vname = apply_aliases(V.gene.name,a);
idx = listmap(vname,TSP.gene.name);
i1 = find(~isnan(idx));
i2 = idx(i1);
V.gene.in_tsp = false(V.ng,1);
V.gene.in_tsp(i1) = true;
V.gene.tsp_n_muts = nan(V.ng,1);
V.gene.tsp_n_muts(i1) = sum(TSP.n_nonsilent(i2,TSP.TOT,:),3);

% correlate to COSMIC data

fprintf('Correlating with COSMIC data\n');
C = load_cosmic_database;
a = load_struct('/xchip/tcga/gbm/analysis/lawrence/tcga/aliases.txt');
C.gene.name = apply_aliases(C.gene.name,a);
vname = apply_aliases(V.gene.name,a);
idx = listmap(vname,C.gene.name);
i1 = find(~isnan(idx));
i2 = idx(i1);
V.gene.in_cosmic = false(V.ng,1);
V.gene.in_cosmic(i1) = true;
V.gene.cosmic_n_muts = nan(V.ng,1);
V.gene.cosmic_n_muts(i1) = C.gene.n_muts(i2);

cd(tmpdir);

end
