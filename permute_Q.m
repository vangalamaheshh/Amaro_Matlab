function QP = permute_Q(Q,D,cyto,by_snp)
  % Randomly permutes segments in matrix Q 
  % Q = segment matrix
  % D = data struct
  % by_snp (default 1): Whether to permute segments by snp or by base
  % pair
    
  if ~exist('by_snp','var') || isempty(by_snp)
    by_snp = 1;
  end
  
  if by_snp
    nsnps = size(D.dat,1);
    
    shift = floor(nsnps*rand(size(Q,1),1));
    QP = Q;
    QP(:,2:3) = mod(QP(:,2:3)+repmat(shift,1,2),nsnps)+1;
    
    cross_chr = find(D.chrn(QP(:,2)) ~= D.chrn(QP(:,3)));
    while ~isempty(cross_chr)
      shift = floor(nsnps*rand(length(cross_chr),1));
      QP(cross_chr,2:3) = mod(QP(cross_chr,2:3)+repmat(shift,1,2),nsnps)+1;
      cross_chr = find(D.chrn(QP(:,2)) ~= D.chrn(QP(:,3)));
    end
  else
    
    chr_lengths = arrayfun(@(x) max([cyto(find([cyto.chrn] == x)).end]), ...
                           1:23);
    clengths = [1 cumsum(chr_lengths)+1];
    QP = Q;
    QP(:,2) = D.pos(QP(:,2));
    QP(:,3) = D.pos(QP(:,3));
    nbps = sum(chr_lengths);
    shift = floor(nbps*rand(size(Q,1),1));
    QP(:,2:3) = mod(QP(:,2:3)+repmat(shift,1,2),nbps)+1;
    cst = arrayfun(@(x) max([find(QP(x,2) >= clengths) 1]),1:size(QP, ...
                                                      1));
    cen = arrayfun(@(x) max([find(QP(x,3) >= clengths) 1]),1:size(QP, ...
                                                      1));
    
    cross_chr = find(cst ~= cen);
    ok = setdiff(1:size(QP,1),cross_chr);
    QP(ok,1) = cst(ok);
    QP(ok,2) = QP(ok,2)-clengths(QP(ok,1))';
    QP(ok,3) = QP(ok,3)-clengths(QP(ok,1))';
    %nls = zeros(1,size(QP,1));
    %for i=1:length(ok)
    %  modi(i,1000)
    %  nls(ok(i)) = length(find_snps(D,QP(ok(i),1),QP(ok(i),2),QP(ok(i), ...
    %                                                    3),0));
    %end
    %cross_chr = union(cross_chr,find(nls <= min_snps));
    
    while ~isempty(cross_chr)
      disp(num2str(length(cross_chr)))
      shift = floor(nbps*rand(length(cross_chr),1));
      QP(cross_chr,2:3) = mod(QP(cross_chr,2:3)+repmat(shift,1,2),nbps)+ ...
          1;
      cst = arrayfun(@(x) max([find(QP(x,2) >= clengths) 1]),cross_chr);
      cen = arrayfun(@(x) max([find(QP(x,3) >= clengths) 1]),cross_chr);
      ok = find(cst == cen);
      QP(cross_chr(ok),1) = cst(ok);
      QP(cross_chr(ok),2) = QP(cross_chr(ok),2)-clengths(QP(cross_chr(ok),1))';
      QP(cross_chr(ok),3) = QP(cross_chr(ok),3)- ...
          clengths(QP(cross_chr(ok),1))';
      cross_chr = cross_chr(find(cst ~= cen));
      %for i=1:length(ok)
      %  modi(i,1000)
      %  nls(ok(i)) = length(find_snps(D,QP(ok(i),1),QP(ok(i),2),QP(ok(i), ...
      %                                                  3),0));
      %end
      %cross_chr = union(cross_chr,find(nls <= min_snps));
    end
  end
    
  
