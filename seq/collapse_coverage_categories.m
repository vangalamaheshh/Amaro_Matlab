function C = collapse_coverage_categories(C,K)

if ~exist('K','var')
  % old way of operating (collapse everything to a single category)
  C.cov = sum(C.cov,3); C.ncat = 1;
  C.totcov = C.cov;
  return;
end

C.totcov = sum(C.cov,3);

nk = slength(K);
c = assign_65x4_to_categ_set(K);
c = max(c,[],3);

newcov = zeros(C.nt,C.ns,nk);
fprintf('Collapsing %d samples: ',C.ns);
for s=1:C.ns, fprintf('%d ',s);
  newcov(:,s,:) = squeeze(C.cov(:,s,:))*c;  % note: matrix multiplication
end, fprintf('\n');

C.cov = newcov;
C.categ = K;
C.ncat = nk;



