function smoothscatter(x,y,P);
% smoothscatter(x,y,P)
% based on:
% SMOOTHHIST2D Plot a smoothed histogram of bivariate data.
%   SMOOTHHIST2D(X,LAMBDA,NBINS) plots a smoothed histogram of the bivariate
%   data in the N-by-2 matrix X.  Rows of X correspond to observations.  The
%   first column of X corresponds to the horizontal axis of the figure, the
%   second to the vertical. LAMBDA is a positive scalar smoothing parameter;
%   higher values lead to more smoothing, values close to zero lead to a plot
%   that is essentially just the raw data.  NBINS is a two-element vector
%   that determines the number of histogram bins in the horizontal and
%   vertical directions.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF) plots outliers in the data as points
%   overlaid on the smoothed histogram.  Outliers are defined as points in
%   regions where the smoothed density is less than (100*CUTOFF)% of the
%   maximum density.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,[],'surf') plots a smoothed histogram as a
%   surface plot.  SMOOTHHIST2D ignores the CUTOFF input in this case, and
%   the surface plot does not include outliers.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF,'image') plots the histogram as an
%   image plot, the default.
%
%   Example:
%       X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%            mvnrnd([0 8], [1 0; 0 5], 2000);
%            mvnrnd([3 5], [5 0; 0 1], 2000)];
%       smoothhist2D(X,5,[100, 100],.05);
%       smoothhist2D(X,5,[100, 100],[],'surf');
%
%   Reference:
%      Eilers, P.H.C. and Goeman, J.J (2004) "Enhancing scaterplots with
%      smoothed densities", Bioinformatics 20(5):623-628.

%   Copyright 2009 The MathWorks, Inc.
%   Revision: 1.0  Date: 2006/12/12
%
%   Requires MATLAB® R14.
%

if nargin < 2, error('requires at least x and y'); end
if nargin > 3, error('max nargin = 3'); end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'lambda',20);
P = impose_default_value(P,'nbins',[100 100]);
P = impose_default_value(P,'outliercutoff',-1);   % don't display outliers
P = impose_default_value(P,'plottype','image');
P = impose_default_value(P,'colormap','blue');
P = impose_default_value(P,'oversaturation_amount',0.1);

% check data
if length(x)~=length(y), error('length(x)~=length(y)'); end
X = [as_column(x) as_column(y)];

% remove bad values
idx = find(isnan(X(:,1)) | isnan(X(:,2)) | isinf(X(:,1)) | isinf(X(:,2)));
X(idx,:) = [];

% prepare
minx = min(X,[],1);
maxx = max(X,[],1);
edges1 = linspace(minx(1), maxx(1), P.nbins(1)+1);
ctrs1 = edges1(1:end-1) + .5*diff(edges1);
edges1 = [-Inf edges1(2:end-1) Inf];
edges2 = linspace(minx(2), maxx(2), P.nbins(2)+1);
ctrs2 = edges2(1:end-1) + .5*diff(edges2);
edges2 = [-Inf edges2(2:end-1) Inf];

[n,p] = size(X);
bin = zeros(n,2);
% Reverse the columns of H to put the first column of X along the
% horizontal axis, the second along the vertical.
[dum,bin(:,2)] = histc(X(:,1),edges1);
[dum,bin(:,1)] = histc(X(:,2),edges2);
H = accumarray(bin,1,P.nbins([2 1])) ./ n;

% Eiler's 1D smooth, twice
G = smooth1D(H,P.lambda);
F = smooth1D(G',P.lambda)';
% % An alternative, using filter2.  However, P.lambda means totally different
% % things in this case: for smooth1D, it is a smoothness penalty parameter,
% % while for filter2D, it is a window halfwidth
% F = filter2D(H,P.lambda);

relF = F./max(F(:));
if P.outliercutoff > 0
    outliers = (relF(P.nbins(2)*(bin(:,2)-1)+bin(:,1)) < P.outliercutoff);
end

% set up colors
nc = 256;
if strcmpi(P.colormap,'hot')
  colormap(hot(nc));
elseif strcmpi(P.colormap,'blue') || strcmp(P.colormap,'bluegreen')
  % (OK)
else
  fprintf('WARNING: unknown colormap "%s": using default "blue"\n',P.colormap);
  P.colormap = 'blue';
end
if strcmpi(P.colormap,'blue')  % default
  scale = [1:-(1/(nc-1)):0]';
  map = [scale scale ones(nc,1)];
%  nc2 = round(nc*(1-P.oversaturation_amount));
%  scale = [1:-(1/(nc2-1)):0]';
%  map = [scale scale ones(nc2,1);zeros(nc-nc2,2) ones(nc-nc2,1)];
  colormap(map);
elseif strcmpi(P.colormap,'bluegreen')  % for print version
  scale = [1:-(1/(nc-1)):0]';
  map = [scale (scale+ones(nc,1))/2 ones(nc,1)];
  colormap(map);
end

% display it
switch P.plottype
case 'surf'
    surf(ctrs1,ctrs2,F,'edgealpha',0);
case 'image'
    imagesc(ctrs1,ctrs2,floor(nc.*relF) + 1,[0 (1-P.oversaturation_amount)]*nc);
    set(gca,'ydir','normal');  % to match orientation of scatter()
    hold on
    % plot the outliers
    if P.outliercutoff > 0
        plot(X(outliers,1),X(outliers,2),'.','MarkerEdgeColor',[.8 .8 .8]);
    end
    hold off
end



%-----------------------------------------------------------------------------
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
% This is a better solution, but takes a bit longer for n and m large
% opts.RECT = true;
% D1 = [diff(E,1); zeros(1,n)];
% D2 = [diff(D1,1); zeros(1,n)];
% Z = linsolve([E; 2.*sqrt(lambda).*D1; lambda.*D2],[Y; zeros(2*m,n)],opts);


%-----------------------------------------------------------------------------
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);
