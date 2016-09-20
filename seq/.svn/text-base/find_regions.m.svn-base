function r = find_regions(T,chr,start,xend)
%
% find_regions(T,chr,start,end)
% find_regions(T,chr,pos)
%
% Given region list T, a structure with required fields:
%    chr
%    start
%    end
% 
% Returns region index list r,
% a list of indices to regions in T that overlap the specified coordinates
%
% Mike Lawrence 2008-05-16
%

require_fields(T, {'chr';'start';'end'});

if ~exist('xend','var')
  xend = start
end

r = find(chr==T.chr & start<T.end & xend>T.start);

