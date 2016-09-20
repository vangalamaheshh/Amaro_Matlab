function M = remove_impossible_mutations(M,method)
if ~exist('method','var'), method = 1; end

if method==1   % slow, with reporting

  % remove mutations in gene/samples/categories with zero coverage
  for g=1:M.ng
    for p=1:M.np
      for c=1:M.TOT-1
        if M.N_cov(g,c,p)==0
          if M.n_silent(g,c,p) | M.n_nonsilent(g,c,p)
            fprintf('Removing impossible mutation(s) in gene %s patient %s\n',...
              M.gene.name{g}, M.patient.name{p});
            M.n_missense(g,c,p)=0;
            M.n_silent(g,c,p)=0;
            M.n_nonsense(g,c,p)=0;
            M.n_indel(g,c,p)=0;
            M.n_splice(g,c,p)=0;
            M.n_nonsilent(g,c,p)=0;
          end
        end
      end
    end
  end
  M.n_missense(:,M.TOT,:) = sum(M.n_missense(:,1:M.TOT-1,:),2);
  M.n_silent(:,M.TOT,:) = sum(M.n_silent(:,1:M.TOT-1,:),2);
  M.n_nonsense(:,M.TOT,:) = sum(M.n_nonsense(:,1:M.TOT-1,:),2);
  M.n_indel(:,M.TOT,:) = sum(M.n_indel(:,1:M.TOT-1,:),2);
  M.n_splice(:,M.TOT,:) = sum(M.n_splice(:,1:M.TOT-1,:),2);
  M.n_nonsilent(:,M.TOT,:) = sum(M.n_nonsilent(:,1:M.TOT-1,:),2);

elseif method==2  % quick

  idx=find(M.N_cov(:)==0);
  M.n_nonsilent(idx)=0;
  M.n_missense(idx)=0;
  M.n_silent(idx)=0;
  M.n_indel(idx)=0;
  M.n_splice(idx)=0;
else
  error('Unknown method\n');
end
