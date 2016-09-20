function [Z,Y]=call_genes(CL,rg,genes,tvals)

for k=1:2 % loop over amp and del
  is_del=k-1;
  
  if is_del
    collapse_method=struct('method','min','find_snps_type',2);
  else
    collapse_method=struct('method','max','find_snps_type',2);
  end

  % X - copy number collapsed to genes
  X=collapse_to_genes(CL,rg,genes,'symb',collapse_method);
  X=rmfield_if_exists(X,{'orig','origidx','gorigidx'});

  if size(tvals{k},2)==1
    tvals{k}=repmat(tvals{k},1,size(CL.dat,2));
  end
  t=[ repmat(-Inf,1,size(CL.dat,2)); tvals{k};  repmat(Inf,1,size(CL.dat,2))];
  
  Y{k}=X;
  if is_del
    X.dat=-X.dat;
  end
  
  Y{k}.dat=nan(size(Y{k}.dat));
  Y{k}.orig=X.dat;
  for i=2:size(t,1)
    Y{k}.dat(X.dat<=repmat(t(i,:),size(Y{k}.dat,1),1) & X.dat>repmat(t(i-1,:),size(Y{k}.dat,1),1))=i-2;
  end
end

Z=Y{2};
Z=rmfield(Z,'orig');

% deletion have a higher precedence to amplifications
% replace the 0 with the amplification scores.
Z.dat=-Z.dat; % make deletions negative
Z.dat(Z.dat==0)=Y{1}.dat(Z.dat==0); % overide 0's with amplfication scores
