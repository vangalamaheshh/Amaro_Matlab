function r=annotate_bed_file(infile,outfile,params)

[r,header_line]=read_bed_file(infile);

switch params.method
 case 'genes'
  rg=params.rg;
  if isfield(params,'all_lesions_file')
    all_lesions=read_table(params.all_lesions_file,{'%s%s%s%s%s%s%s%s%s',9,'%f'},char(9),1);
  else
    all_lesions=[];
  end

  if isfield(params,'region_annotation_db')
    region_annotation_db=read_table(params.region_annotation_db,'%s%f%f%s%s',char(9),1);
%    region_annotation_db.dat{1}=chromosome2num(region_annotation_db.dat{1});    
  else
    region_annotation_db=[];
  end
  
  for i=1:slength(r)
    disp(i);
    genes='';
    ampdel=0;
    if isfield(params,'ampdel')
      if ~isempty(grep(params.ampdel{1},r.name{i},1))
        ampdel=1;
      elseif ~isempty(grep(params.ampdel{2},r.name{i},1))
        ampdel=2;
      end
    end
    if ampdel>1 % if 2 == del
      if isfield(params,'gene_gistic') && params.gene_gistic==1
       [genes_in_region,closest_to_region]=genes_at(rg,chromosome2num(r.chr{i}),r.start(i),r.end(i),1,0);
      else
       [genes_in_region,closest_to_region]=genes_at(rg,chromosome2num(r.chr{i}),r.start(i),r.end(i),1,1);
      end
    else
      [genes_in_region,closest_to_region]=genes_at(rg,chromosome2num(r.chr{i}),r.start(i),r.end(i),1,1);
    end      
    %%% unique the above by locuslink i
    if ~isempty(genes_in_region)
      lid_in_region = cat(1,rg(genes_in_region).locus_id);
      [un_lid_reg,un_lid_reg_i1,un_lid_reg_i2] = unique(lid_in_region);
      genes_in_region=genes_in_region(un_lid_reg_i1);
    end    
    
    if length(genes_in_region)==0
      % negative indicates closest_to_region
      gene_ids=-closest_to_region;
      gene_symbs=rg(closest_to_region).symb;
      
      if isfield(params,'gene_lists') && isfield(params,'markers')
        gene_annot_symbs=annotate_gene_list(gene_symbs,params.gene_lists,params.markers);
        gene_annot_symbs=gene_annot_symbs{1};
      end
      
      r.n_genes(i)=0;
      if isfield(params,'E')
        Egidx=grep(['^' gene_symbs '$'],params.E.gacc,1);
        Cgidx=grep(['^' gene_symbs '$'],params.C.gacc,1);
        if ~isempty(Cgidx)
          tmp=histc(params.C.dat(Cgidx(1),:),params.hist_edges);
          if params.hist_trim_last
            tmp(end)=[];
          end
          r.copy_counts(i,:)=tmp;
        else
          r.copy_counts(i,:)=nan(1,5);
        end
        if ~isempty(Egidx) && ~isempty(Cgidx)
          c=corrcoef(params.E.dat(Egidx(1),:),params.C.dat(Cgidx(1),:));
          cstr=num2str(c(1,2));
        else
          cstr='NA';
        end
        genes=sprintf('[%s(%s)]',gene_annot_symbs,cstr);    
      else
        genes=sprintf('[%s]',gene_annot_symbs);    
      end
    else
      r.n_genes(i)=length(genes_in_region);
      gene={};
      gene_ids=genes_in_region;
      
      gene_symbs={rg(genes_in_region).symb};
      
      if isfield(params,'gene_lists') && isfield(params,'markers')
        gene_annot_symbs=annotate_gene_list(gene_symbs,params.gene_lists,params.markers);
      end
      
      cs=[];
      for gi=1:length(genes_in_region)
%        disp(genes_in_region(gi));
        if isfield(params,'E')
          Egidx=grep(['^' gene_symbs{gi} '$'],params.E.gacc,1);
          Cgidx=grep(['^' gene_symbs{gi} '$'],params.C.gacc,1);
          if ~isempty(Egidx) && ~isempty(Cgidx)
            c=corrcoef(params.E.dat(Egidx(1),:),params.C.dat(Cgidx(1),:));
            cs(gi)=c(1,2);
            gene{gi}=sprintf('%s(%4.2f),',gene_annot_symbs{gi},cs(gi));
          else
            cs(gi)=NaN;
            gene{gi}=sprintf('%s(NA),',gene_annot_symbs{gi});
          end
        else
          gene{gi}=sprintf('%s,',gene_annot_symbs{gi});
        end
      end
      

      % sort by correlation
      if ~isempty(cs)
        cstmp=cs;
        cstmp(isnan(cs))=-100;
        [css,csi]=sort(cstmp,'descend');
      else
        csi=1:length(genes_in_region)
      end
    
      gene_ids=gene_ids(csi);
      
      % build string
      genes=[];
      for gi=1:length(genes_in_region)
        genes=[genes gene{csi(gi)}];
      end
      
      % calc copy_counts
      Cgidx=grep(['^' rg(genes_in_region(csi(1))).symb '$'],params.Call.gacc,1);
      if ~isempty(Cgidx)
        tmp=histc(params.Call.dat(Cgidx(1),:),params.hist_edges);
        if params.hist_trim_last
          tmp(end)=[];
        end
        r.copy_counts(i,:)=tmp;
      else
        r.copy_counts(i,:)=nan(1,5);
      end
      
    end
    r.genes{i}=genes;
    r.gene_ids{i}=gene_ids;
   
    if ~isempty(all_lesions)
      idx=grep(sprintf('%s:%d-%d',r.chr{i},r.start(i),r.end(i)),all_lesions.dat{3},1);
      if isempty(idx)
        disp(['Could not find ' num2str(i)]);
      else
        st=sprintf('%d\t',idx(1));
        for ai=1:8
          st=[st sprintf('%s\t',all_lesions.dat{ai}{idx(1)})];
        end
        st(end)=[]; % remove last tab
        
        r.all_lesions_str{i}=st;
      end
    end
    
    if ~isempty(region_annotation_db)
      idx=intersect(strmatch(r.chr{i},region_annotation_db.dat{1}),...
                    find(region_annotation_db.dat{3}>= r.start(i) & ...
                         region_annotation_db.dat{2}<= r.end(i)));
      if isempty(idx)
        r.region_annotation_db{i}=sprintf('\t');
      else
        r.region_annotation_db{i}=sprintf('%s\t%s',region_annotation_db.dat{4}{idx(1)}, ...
            region_annotation_db.dat{5}{idx(1)});
      end     
    end
    
  end
end


r.genes=r.genes';
r.gene_ids=r.gene_ids';
r.n_genes=r.n_genes';
r.all_lesions_str=r.all_lesions_str';
r.region_annotation_db=r.region_annotation_db';

f=fopen(outfile,'w');
fprintf(f,'%s\n',header_line);
for i=1:slength(r)
  fprintf(f,'%s\t%d\t%d\t%s',r.chr{i},r.start(i),r.end(i),r.name{i});
  if isfield(r,'all_lesions_str')
    fprintf(f,'\t%s',r.all_lesions_str{i});
  end
  if isfield(r,'region_annotation_db')
    fprintf(f,'\t%s',r.region_annotation_db{i});
  end
  if isfield(r,'other')
    for j=1:size(r.other,2)
      fprintf(f,'\t%s',r.other{i,j});
    end
  end
  fprintf(f,'\t%d\t%d\t%d\t%d\t%d',r.copy_counts(i,:));
  fprintf(f,'\t%d',r.n_genes(i));
  fprintf(f,'\t%s',r.genes{i});
  fprintf(f,'\n');
end
fclose(f);
