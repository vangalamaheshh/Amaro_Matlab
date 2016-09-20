function Tc = flank_regions(T, bp)
%
% flank_regions(T, bp)
%
% Given region list T, a structure with required fields:
%    gene
%    chr
%    strand
%    start
%    end
% 
% Returns region list Tc
% in which each regions has been extended by "bp" bases in each direction.
%
% Mike Lawrence 2008-05-15
%

require_fields(T, {'gene';'chr';'strand';'start';'end'});

Tc = T;
for i=1:length(Tc.gene)
  Tc.start(i) = Tc.start(i) - bp;
  Tc.end(i) = Tc.end(i) + bp;
end

