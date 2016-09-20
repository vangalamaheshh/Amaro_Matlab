function write_gene_table(fname,G,cyto,rg)
%WRITE_GENE_TABLE save gene-level copy number data to file
%
%   WRITE_GENE_TABLE(FNAME,G,CYTO,RG)
%
% Outputs G (copy number reduced to gene level) to the file FNAME, 
% using CYTO and RG to add locus ID and cytoband annotation columns.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

  nsamples = size(G.dat,2);
  ngenes = size(G.dat,1);
  f = fopen(fname,'w');

  % print header
  fprintf(f,'%s\t%s\t%s','Gene Symbol','Locus ID','Cytoband');
  for j=1:nsamples
    fprintf(f,'\t%s',char(G.sdesc(j,:)));
  end
  fprintf(f,'\n');
  
  fmt = ['%s\t%d\t%s' repmat('\t%1.3f',1,nsamples) '\n'];
  for j=1:ngenes
    modi(j,1000)
%!    if G(j).chrn < 23 %!!!HUMCHR
      band = find([cyto.chrn] == G.chrn(j) & [cyto.start] <= G.gstart(j),1,'last');
      locus_id = [];
      if isfield(rg,'locus_id')
          % (if multiple locus IDs, use the first)
          locus_id = rg(G.rgindex{j}(1)).locus_id;
      end
      fprintf(f,fmt,G.gdesc{j},locus_id,cyto(band).name,G.dat(j,:)); 
%!    end
  end  
  fclose(f);

