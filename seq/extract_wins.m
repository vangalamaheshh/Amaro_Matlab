function x = extract_wins(g,windowsize)
if ~exist('windowsize','var'), windowsize=100;end

clen = load_chrlen;
nwins = round(1.1*clen/windowsize);

x = cell(24,1);
for c=1:24
  x{c} = g.getContents(c,1,nwins(c));
end
x = cat(1,x{:});
x(x==-128)=0;  % empty positions




