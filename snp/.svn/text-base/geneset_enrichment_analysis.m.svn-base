function [M chis ps F num_regs genes_in_list] = geneset_enrichment_analysis(outfile,geneset_desc_file,D,regs,cyto,rg,ampdel,gs_size_range,max_genes_per_reg,max_n_per_reg,rm_provisional)
  
  % geneset descriptor file should containing following files
    % Ds: n x 2 cell array; first column has names of geneset, second
    % column has geneset description
    % gene_lists: genes contained in each geneset
    if ~exist('max_n_per_reg','var') 
      max_n_per_reg = [];
    end
    
    if ~exist('rm_provisional','var') || isempty(rm_provisional)
      rm_provisional = 0;
    else
        if rm_provisional
          status = {rg.status};  
          predicted_idx = strmatch('Predicted',status);
          inferred_idx = strmatch('Inferred',status);
          provisional_idx = strmatch('Provisional',status);
          rm_idx = cat(1,predicted_idx,inferred_idx,provisional_idx);
          keep_idx = setdiff(1:size(rg,2),rm_idx);
          rg = rg(keep_idx);
        end
    end
        
    %% Load geneset descriptor file files
    verbose('Loading files...',20)
    load(geneset_desc_file);
    
    num_genes = cellfun(@length,gene_lists);
    
    %load /kinome/snp/snp01/Craig/GCM/gene_lists/go_defs.mat
    
    % Add genes to regs
    if ~isfield(regs{1},'genes')
      regs = add_genes(D,regs,cyto,rg);
    end
    
    % Select gos to analyze
    
    keep_gs = find(num_genes >= gs_size_range(1) & num_genes <= ...
                   gs_size_range(2));
    
    % Do analysis
    
    M = zeros(2,3,length(keep_gs));
    chis = zeros(1,length(keep_gs));
    ps = zeros(1,length(keep_gs));
    F = zeros(2,3,length(keep_gs));
    num_regs = cell(1,length(keep_gs));
    genes_in_list = cell(1,length(keep_gs));
    
    for j=1:length(keep_gs)
      if mod(j,10) == 0
        verbose(num2str(j),20)
      end
      [M(:,:,j) chis(j) ps(j) F(:,:,j) num_regs{j} genes_in_list{j}] = gen_geneset_contingency_table(regs,rg,[],max_genes_per_reg,ampdel,gene_lists{keep_gs(j)},[],[],max_n_per_reg);
    end
    
    qs = calc_fdr_value(ps);
    
    [S si] = sort(chis,'descend');
    
    
    % Write output
    
    f = fopen(outfile,'w');
    
    fprintf(f,'Rank');
    fprintf(f,'\t');
    fprintf(f,'GO term');
    fprintf(f,'\t');
    fprintf(f,'Description');
    fprintf(f,'\t');
    fprintf(f,'Chi-value');
    fprintf(f,'\t');
    fprintf(f,'Chi-p');
    fprintf(f,'\t');
    fprintf(f,'Chi-q');
    fprintf(f,'\t');
    fprintf(f,'Num Genes in GO term');
    ampst = {'Amp','Del'};
    for k=1:2
      if ampdel(k)
      fprintf(f,'\t');
      fprintf(f,['Number of ' ampst{k} ' regs containing gene in GO category.']);
      fprintf(f,'\t');
      fprintf(f,[ampst{k} ' genes in category']);
      end
    end
    fprintf(f,'\n');
    
    
    for j=1:length(si)
      if mod(j,100) == 0
        verbose(num2str(j),20)
      end
      cur_go = Ds(keep_gs(si(j)),1);
      cur_desc = char(Ds(strmatch(cur_go,Ds),2));
      fprintf(f,num2str(j));
      fprintf(f,'\t');
      fprintf(f,char(cur_go));
      fprintf(f,'\t');
      fprintf(f,cur_desc);
      fprintf(f,'\t');
      fprintf(f,'%3.4f',chis(si(j)));
      fprintf(f,'\t');
      fprintf(f,'%0.3g',ps(si(j)));
      fprintf(f,'\t');
      fprintf(f,'%0.3g',qs(si(j)));
      fprintf(f,'\t');
      fprintf(f,num2str(num_genes(keep_gs(si(j)))));
      for k=1:2
        if ampdel(k)
          fprintf(f,'\t');
          cur_regs = num_regs{si(j)};
          fprintf(f,num2str(cur_regs(k)));
          fprintf(f,'\t');
          reggenes = genes_in_list{si(j)};
          cur_genes = reggenes{k};
          str = '';
          for q = 1:length(cur_genes)
            str = [str ' ' cur_genes{q} ';'];
          end
          fprintf(f,str(1:end-1));
        end
      end
      fprintf(f,'\n');
    end
    
    fclose(f);
