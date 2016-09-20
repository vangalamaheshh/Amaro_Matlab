function [z c] = factor2lego(h,catname)
% [z c] = factor2lego(h,catname)
%
% Helper function for LEGO-plot visualization of NMF factors
%
% GIVEN "h", a set of factor definitions from NMF as performed in mutation_spectra_pca_plot.m
%          each row = a factor
%          columns = the 96 categories of mutations (after strand collapse)
% and "catname"
%         = tells what order the 96 categories are in.  (default order used if not provided)
%
% REORDERS "h" to yield:
%       z = a set of 8x12 heights for LEGO plot
%
% RETURNS also:
%       c = a set of 8x12 colors for LEGO plot
%
% Note: if "h" has multiple rows, then "z" will have multiple pages
%       (but c will not be provided in duplicate)
%
% Mike Lawrence 2012-05-22

if ~exist('catname','var')
  catname = get_default_catname();
end

if size(h,2)~=96, error('h should have 96 columns'); end

classes = {'C->G','C->A','C->T','A->T','A->C','A->G'};
bases = 'TCAG';

z = zeros(8,12,size(h,1));
cx=1; cy=1;
for ci=1:6
  for left=1:4
    for right=1:4
      name = [classes{ci}(1) ' in ' bases(left) '_' bases(right) ' ' classes{ci}(2:end)];
      cidx = find(strcmp(catname,name));
      if length(cidx)~=1, error('Problem with catname'); end
      y = cy+left-1;
      x = cx+right-1;
      z(y,x,:) = h(:,cidx);
  end,end
  cx=cx+4;
  if cx>9, cx=1;cy=cy+4; end
end

if nargout>=2
  c = get_LEGO_colors();
end























end


function x = get_default_catname()
  x = {...
    'C in A_A ->A'
    'C in A_C ->A'
    'C in A_G ->A'
    'C in A_T ->A'
    'C in C_A ->A'
    'C in C_C ->A'
    'C in C_G ->A'
    'C in C_T ->A'
    'C in G_A ->A'
    'C in G_C ->A'
    'C in G_G ->A'
    'C in G_T ->A'
    'C in T_A ->A'
    'C in T_C ->A'
    'C in T_G ->A'
    'C in T_T ->A'
    'A in A_A ->C'
    'A in A_C ->C'
    'A in A_G ->C'
    'A in A_T ->C'
    'A in C_A ->C'
    'A in C_C ->C'
    'A in C_G ->C'
    'A in C_T ->C'
    'A in G_A ->C'
    'A in G_C ->C'
    'A in G_G ->C'
    'A in G_T ->C'
    'A in T_A ->C'
    'A in T_C ->C'
    'A in T_G ->C'
    'A in T_T ->C'
    'A in A_A ->G'
    'A in A_C ->G'
    'A in A_G ->G'
    'A in A_T ->G'
    'A in C_A ->G'
    'A in C_C ->G'
    'A in C_G ->G'
    'A in C_T ->G'
    'A in G_A ->G'
    'A in G_C ->G'
    'A in G_G ->G'
    'A in G_T ->G'
    'A in T_A ->G'
    'A in T_C ->G'
    'A in T_G ->G'
    'A in T_T ->G'
    'C in A_A ->G'
    'C in A_C ->G'
    'C in A_G ->G'
    'C in A_T ->G'
    'C in C_A ->G'
    'C in C_C ->G'
    'C in C_G ->G'
    'C in C_T ->G'
    'C in G_A ->G'
    'C in G_C ->G'
    'C in G_G ->G'
    'C in G_T ->G'
    'C in T_A ->G'
    'C in T_C ->G'
    'C in T_G ->G'
    'C in T_T ->G'
    'A in A_A ->T'
    'A in A_C ->T'
    'A in A_G ->T'
    'A in A_T ->T'
    'A in C_A ->T'
    'A in C_C ->T'
    'A in C_G ->T'
    'A in C_T ->T'
    'A in G_A ->T'
    'A in G_C ->T'
    'A in G_G ->T'
    'A in G_T ->T'
    'A in T_A ->T'
    'A in T_C ->T'
    'A in T_G ->T'
    'A in T_T ->T'
    'C in A_A ->T'
    'C in A_C ->T'
    'C in A_G ->T'
    'C in A_T ->T'
    'C in C_A ->T'
    'C in C_C ->T'
    'C in C_G ->T'
    'C in C_T ->T'
    'C in G_A ->T'
    'C in G_C ->T'
    'C in G_G ->T'
    'C in G_T ->T'
    'C in T_A ->T'
    'C in T_C ->T'
    'C in T_G ->T'
    'C in T_T ->T'
  };
end
