function [GrX,CrX_genes]=outlier_expression(Gr,Cr,chr,rng,ts_x1,ts_x2,amp_or_del)

% keep genes in range
grange=cat(1,Gr.grange{:});
with_grange=find(~any(isnan(grange),2));

genes_in_range=find(Gr.chrn==chr & grange(:,1)<=rng(2) & grange(:,2)>=rng(1));
GrX=reorder_D_rows(Gr,genes_in_range);
if ~isempty(genes_in_range)

  % sort by position
  mid=mean(cat(1,GrX.grange{:}),2);
  [sm,smi]=sort(mid);
  GrX=reorder_D_rows(GrX,smi);
  
  % get copy number per genes
  length(smi)
  pos=chrpos2genomic_location(GrX.chrn,cat(1,GrX.grange{:}));
  [CrX_genes,v_genes]=collapse_to_genes(Cr,[],pos);
  
  c=CrX_genes.dat;
  e_norm=GrX.dat;
  e_low=e_norm;
  
  e_norm(c>ts_x1(1) | c<-ts_x1(2))=NaN; 
  switch amp_or_del
   case 1
    e_low(c<=ts_x1(1))=NaN; % ignore non-amp or deleted samples
    disp('here');
    %  e_low(c>ts_x2(1))=NaN; % ignore high level of amplification
   case 2
    e_low(c>=-ts_x1(2))=NaN; % ignore non-amp or amplified samples
                             %  e_low(c<-ts_x2(2))=NaN; % ignore HD
  end
  mmed=nanmedian(e_norm,2);
  if (0)
    mmad=nanmedian( abs([e_low-repmat(nanmedian(e_low,2),1,size(e_low,2)) e_norm-repmat(nanmedian(e_norm,2),1,size(e_norm,2)) ]),2);
%    mmad=nanmedian( abs([ e_norm-repmat(nanmedian(e_norm,2),1,size(e_norm,2)) ]),2);
  else
    mmad=mad(e_norm,1,2);
  end
  % mmad(mmad<0.1)=0.1;
  % keyboard
  
  switch amp_or_del
   case 1
    p=prctile(e_low,90,2);
    d=p-mmed;
    z=d./mmad; 
    GrX=add_D_sup(GrX,'OSPRC','Outlier Score - prctile',p','rows');
    GrX=add_D_sup(GrX,'OSDEL','Outlier Score - delta',d','rows');
   case 2
    p=prctile(e_low,10,2);
    d=-(p-mmed);
    z=d./mmad;
    GrX=add_D_sup(GrX,'OSPRC','Outlier Score - prctile',p','rows');
    GrX=add_D_sup(GrX,'OSDEL','Outlier Score - delta',d','rows');
  end
  
  GrX=add_D_sup(GrX,'OSMED','Outlier Score - median',mmed','rows');
  GrX=add_D_sup(GrX,'OSMAD','Outlier Score - mad',mmad','rows');
  GrX=add_D_sup(GrX,'OS','Outlier Score',z','rows');
  
  [sz,szi]=sort(-z);
  GrX=reorder_D_rows(GrX,szi);
end
  CrX_genes=[];
end



