function qual = dRanger_calculate_quality(X,sparams)
% dRanger_calculate_quality(X,sparams)
%
% computes "quality" (0-1) of each rearrangement
%   from the measured locus characteristics of end1 and end2

% input fields

flds = {'tumreads','fmapqzT1','fmapqzN1','fmapqzT2','fmapqzN2','nuwpT1',...
  'nuwpN1','nuwpT2','nuwpN2','zstdev1','zstdev2'};
demand_fields(X,flds);
X = make_numeric(X,flds);

% if normpanel was used, find out what the maximum # samples were
normpanelsize = 0;
if isfield(X,'normpaneldetails')
  z=[]; z.x = grep('^(\d+)_samples:0$',unique(X.normpaneldetails));
  z = parse_in(z,'x','^(\d+)_samples:0$','num',1);
  normpanelsize = max(z.num);
end
 
% sigmoidal curve with three parameters:
%    M = min value for penalty to be imposed
%    K = inflection point
%    Z = sigmoidicity

if ~exist('sparams','var')
  
  if normpanelsize<10
    %          M    K    Z
    sparams = [1    1    3;  ...       % z-score of stdev
               0.25 0.25 2;  ...       % frac mapq=0
               4    4    3   ];        % # unique weirdpairs

  else  % with large enough panel of normals, disable the mapqz penalty
    %          M    K    Z
    sparams = [1    1    3;  ...       % z-score of stdev
               0.95 0.98 2;  ...       % frac mapq=0
               4    4    3   ];        % # unique weirdpairs
  end
end

testtype = [1 1 2 2 2 2 3 3 3 3];
t.M = sparams(testtype,1)';
t.K = sparams(testtype,2)';
t.Z = sparams(testtype,3)';

% for poorly supported rearrangements (<4 pairs), increase stringency for low stdevs

penalty = 0 + (1.25*(X.tumreads==3)) + (1.5*(X.tumreads==2));
zs1 = X.zstdev1; neg = find(zs1<0); zs1(neg) = zs1(neg).*penalty(neg); zs1 = abs(zs1);
zs2 = X.zstdev2; neg = find(zs2<0); zs2(neg) = zs2(neg).*penalty(neg); zs2 = abs(zs2);
m1 = [zs1 zs2];

% for extremely well-supported rearrangements, relax stringency of fmapqz/nuwp filtering
relax = 1-0.9*(min(500,X.tumreads)/500);
relax(X.tumreads>500) = 0.1;
relax(X.tumreads<20) = 1;
m2 = [X.fmapqzT1 X.fmapqzN1 X.fmapqzT2 X.fmapqzN2 X.nuwpT1 X.nuwpN1 X.nuwpT2 X.nuwpN2];
m2 = bsxfun(@times,m2,relax);

% compute total quality
metrics = [m1 m2];
qual = nan(slength(X),1);
for i=1:slength(X)
  tmp = metrics(i,:)-t.M;
  tmp = max(0,tmp);
  tmp = tmp ./ t.K;
  tmp = tmp .^ t.Z;
  tmp = 1 - (tmp ./ (1+tmp));
  tmp = prod(tmp);
  qual(i) = tmp;
end

