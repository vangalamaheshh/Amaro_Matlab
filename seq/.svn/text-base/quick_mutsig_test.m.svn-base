function G = quick_mutsig_test(G,P);

ng = slength(G);

if ~exist('P','var'), P=[];end
P = impose_default_value(P,'genes_to_analyze',1:ng);
P = impose_default_value(P,'report_top_n_mutsig_genes',min(25,length(P.genes_to_analyze)));

if iscellstr(P.genes_to_analyze)
  demand_field(G,'name');
  P.genes_to_analyze = listmap(P.genes_to_analyze,G.name);
end
analyze_this_gene = ismember(1:ng,P.genes_to_analyze);

globalrate_non = sum(G.nnon)/sum(G.Nnon);

if ~isfield(G,'Fobs'), G.Fobs = (G.nnon./G.Nnon)/globalrate_non; end

if isfield(G,'Ffit_distrib')

  P=impose_default_value(P,'newmethod1_Frange_to_test','*required*');
  Frange = P.newmethod1_Frange_to_test;

  G.Ffit = nan(slength(G),1);
  G.p = nan(slength(G),1);
  for gi=1:slength(G), modi(gi,100);
    if analyze_this_gene(gi)
      Fdist = G.Ffit_distrib{gi};
      pmutsig = 1-binocdf(G.nnon(gi)-1,G.Nnon(gi),globalrate_non*Frange);
      pmutsig_weighted = sum(pmutsig.*Fdist);
      G.p(gi) = pmutsig_weighted;
      [tmp idx] = max(Fdist); G.Ffit(gi) = Frange(idx);

      if 0&&strcmp(G.name{gi},'KEAP1')
        y = [Fdist pmutsig];
        y = bsxfun(@rdivide,y,max(y)*1.1);
        plot(Frange,y); title(G.name{gi},'fontsize',20);set(gca,'ytick',[]);
        pr(Frange,Fdist,pmutsig,pmutsig.*Fdist,cumsum(pmutsig.*Fdist));
      end

    end
  end

elseif isfield(G,'Ffit')

  localrate = globalrate_non * G.Ffit;
  G.p = 1-binocdf(G.nnon-1,max(G.nnon,G.Nnon),localrate); G.p(G.Nnon==0) = 1;

end

G.q = calc_fdr_value(G.p);

G.ratio = G.Fobs./G.Ffit;

[tmp ord] = sort_struct(G,{'p','Fobs'},[1 -1]); [tmp G.rank] = sort(ord);
flds = {'rank','name','Nnon','nsil','nnon','Fsil','Fobs','Nfit','nfit','Ffit','ratio','p','q'};

pr(keep_fields_that_exist(G,flds),ord(1:P.report_top_n_mutsig_genes));

if nargout==0, clear G; end
