function B=read_by_gene_file(fname,rg)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

tbl=read_table(fname,{'%s%s%s',3,'%f'},char(9),1,'commentstyle','shell');
B.sdesc=tbl.headers{1}(4:end);
B.dat=cat(2,tbl.dat{4:end});
B.gacc=tbl.dat{1};
B.gdesc=tbl.dat{2};
B.cyto=tbl.dat{3};

if exist('rg','var')
  if ischar(rg)
    rgfname=rg;
    load(rgfname);
  end
  [Mt,m1,m2]=match_string_sets_hash({rg.symb},B.gacc);
  [um2,umi,umj]=unique(m2);
  um1=m1(umi);
  B.chr={rg(um1).chr};
  B.chrn=chromosome2num(B.chr);
  B.pos=floor(0.5*(cat(1,rg(um1).start)+cat(1,rg(um1).end)));
end
