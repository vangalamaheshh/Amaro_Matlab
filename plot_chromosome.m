function plot_chromosome(ch,cyto,segments,priority,orientation,disp_p,in_pos)


if exist('in_pos','var') && ~isempty(in_pos)
  gr=make_subplotgrid(disp_p.x.sizes,disp_p.y.sizes,disp_p.x.gaps, ...
                      disp_p.y.gaps,disp_p.x.border,disp_p.y.border,in_pos);
else
  gr=make_subplotgrid(disp_p.x.sizes,disp_p.y.sizes,disp_p.x.gaps, ...
                      disp_p.y.gaps,disp_p.x.border,disp_p.y.border);  
end

