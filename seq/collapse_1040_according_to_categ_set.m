function c = collapse_1040_according_to_categ_set(C)
% c = collapse_1040_according_to_categ_set(C)
%
% Given a set of categories as struct C with the following fields:
%   left   = subset of 'ACGT', representing 5' base
%   right  = subset of 'ACGT', representing 3' base
%   from   = subset of 'AC', representing mutated base (after strand collapse)
%
% Maps the full set of 1040 territory categories onto this reduced set of categories
%   returns zero for catgories that don't map (e.g. "exon:any N")
%
% Mike Lawrence 2010-01-27

require_fields(C,{'left','right','from'});

X = load_struct('/xchip/tcga_scratch/lawrence/db/allcateg/categs.txt','%f%s');
require_fields(X,{'num','name'});
X = sort_struct(X,'num');  % (already sorted, but doesn't hurt to make sure)

c = zeros(slength(X),1);

x = parse(X.name,'\S+:\S+:(.) in (.)_(.)'),{'from','left','right'});


%%%%%%%%%%%%%%% (NOT FINISHED)


% (code from melanoma project)

X = load_struct('/xchip/tcga_scratch/lawrence/db/allcateg/categs.txt');
X.categ5 = zeros(slength(X),1);
X.categ5(grep('(A|T) in',X.name,1)) = 5;
X.categ5(grep('C in (C|T)_(C|G)',X.name,1)) = 1;
X.categ5(grep('G in (C|G)_(A|G)',X.name,1)) = 1;
X.categ5(grep('C in (C|T)_(A|T)',X.name,1)) = 2;
X.categ5(grep('G in (A|T)_(A|G)',X.name,1)) = 2;
X.categ5(grep('C in (A|G)_(A|T)',X.name,1)) = 3;
X.categ5(grep('G in (A|T)_(C|T)',X.name,1)) = 3;
X.categ5(grep('C in (A|G)_(C|G)',X.name,1)) = 4;
X.categ5(grep('G in (C|G)_(C|T)',X.name,1)) = 4;
save_textfile(sprintf('%.0f\n',X.categ5),'/xchip/tcga_scratch/lawrence/mel/analysis/20091205/terr_collapse.txt');

