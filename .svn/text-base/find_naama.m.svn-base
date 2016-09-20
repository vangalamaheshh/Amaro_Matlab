function [sc,sr]=find_naama(dat,rowcore,row_std,col_std,R,C)
verbose=0;
if nargin == 4
  R=dna_norm(dat);
  C=dna_norm(dat')';
end
%num_std=3;
sr=[]; sc=[];
coreC=sum(C(rowcore,:),1);
sru=find(coreC > col_std * sqrt(length(rowcore)) );
%srd=find(coreC <  -col_std * sqrt(length(rowcore)) );
%if ( length(srd) > length(sru) )
%  sr=srd;
%else
  sr=sru;
%end
%sr=find(abs(coreC) > num_std * sqrt(length(rowcore)) );
if(verbose)
  figure(1)
  subplot(2,1,1);
  comp_dist({coreC(sr)',coreC(setdiff(1:length(coreC),sr))'},50);
end
if ~isempty(sr) 
  srR=sum(R(:,sr).*repmat(sign(coreC(:,sr)),size(dat,1),1),2);
  sc=find(srR > row_std * sqrt(length(sr)) );

  if(verbose)
    subplot(2,1,2);
    comp_dist({srR(sc)',srR(setdiff(1:length(srR),sc))},50);
  end
%  [mean(coreC) std(coreC); ...
%   mean(coreC(:,sr)) std(coreC(:,sr)); ...
%   mean(srR) std(srR); ...
%   mean(srR(sc)) std(srR(sc)); ...
%  ]
else
  warning('Did not find any rows');
end
% pause




