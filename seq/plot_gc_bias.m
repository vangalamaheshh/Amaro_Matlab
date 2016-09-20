function plot_gc_bias(gc,cov,gcwin,covcap)

if ~exist('gcwin','var'), fprintf('Assuming gcwin=200\n'); gcwin=200; end
if ~exist('covcap','var'), fprintf('Using covcap=10000\n'); covcap=10000; end

n = min(length(gc),length(cov));
gc = gc(1:n);
cov = cov(1:n);

H = zeros(covcap+1,gcwin+1);

for i=1:n
  if ~mod(i,100000), fprintf('%d/%d ',i,n); end
  y = min(covcap,cov(i));
  if y<0, error('coverage should be nonnegative'); end
  x = min(gcwin,gc(i));
  if x==-1, continue; end   % N regions of genome
  H(y+1,x+1) = H(y+1,x+1) + 1;
end, fprintf('\n');

imagesc(H); colorbar;
ylabel('coverage'); xlabel('GC');

fprintf('Type "return" to exit function\n');
keyboard
