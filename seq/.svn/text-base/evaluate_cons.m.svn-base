function cons = evaluate_cons(c,st,en)
%
% evaluate_cons(c,st,en)
%
% c is the uint8 array loaded using load_cons(chr)
% st:en is the region across which to take the average conservation
%
% returns an array of conservation values
%    ranging from 0.000 to 1.000
%    values absent from data files are given as NaN
%    values out of range are given as NaN
%

if st<1, error('st must be at least 1'); end

cons = nan(en-st+1,1);
en = min(en,length(c));
cons(1:en-st+1) = double(c(st:en));
cons(cons==255)=nan;
cons = cons/250;
