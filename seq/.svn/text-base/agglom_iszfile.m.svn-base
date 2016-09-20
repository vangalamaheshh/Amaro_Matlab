function agglom_iszfile(infile,outfile,P)
% given an isz file
%     <sz> <readgroup1> <readgroup2> <readgroup3> ... <reagroupn>
%
% identifies readgroups that are similar enough to collapse, and sums all their counts
% identifies readgroups that have too few counts to identify, and zeros all their counts
%
% writes an isz file with the same number of columns and rows
%
% Mike Lawrence 2009-12-28

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'similarity_cutoff',0.0007);
P = impose_default_value(P,'min_lane_nreads',100000);

tmp = load_matrix(infile);
sz = tmp(:,1);
nsz = size(tmp,1);
nln = size(tmp,2)-1;
I = tmp(:,2:end);
lanetot = sum(I,1);
In = bsxfun(@rdivide,I,lanetot);

d = dist(In',[],'euclidean');
dn = bsxfun(@rdivide,d,sqrt(nsz));

id = sparse(dn<P.similarity_cutoff);
[S C] = graphconncomp(id,'Directed',false);
C(lanetot<P.min_lane_nreads)=0;

O = zeros(nsz,nln);
[u ui uj] = unique(C);
for i=1:length(u)
  if u(i)==0, continue; end   % lanes with too few counts are zeroed out
  idx = find(uj==i);
  O(:,idx) = repmat(sum(I(:,idx),2),1,length(idx));
end

save_matrix([sz O],outfile);

return

% visualize the transform
olanetot = sum(O,1);
On = bsxfun(@rdivide,O,olanetot);
figure(1),clf
subplot(2,1,1), plot(In), xlim([0 500]);
subplot(2,1,2), plot(On), xlim([0 500]);
keyboard
