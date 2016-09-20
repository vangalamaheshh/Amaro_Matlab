function fet_w=weigh_features(D,method,prior)

[ord,pv,st]=gsea_get_order(D.dat,D.supdat(1,:),method.type);
fet_w=-log10(pv);

if isfield(method,'select')
  switch method.select
   case 'topk'
    [ss,si]=sort(fet_w);
    fet_w(si(1:(end-method.k)))=0;
  end
end

