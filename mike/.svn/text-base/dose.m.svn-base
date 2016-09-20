% dose.m

% Mike Lawrence 2008-01-05

% Given matched matrices of copy-number and expression data,
% computes correlation between expression and copy number for each gene,
% and plots the distribution of correlations as a histogram.
% Also computes a matrix of "trans effect" residuals after
% linear regression to subtract "cis effect" of gene dosage.
% Generates bitmap representations of input and output matrices.
% Outputs regression parameters and residuals as tablular files.

%%% Parameters

% input files

cfilename = 'Copy_number_11776_matched_to_EX.20080102.txt';
efilename = 'Expression_11776_matched_to_CN.20080102.txt';

ns = 80;
ng = 11776;

% output files

pfilename = 'Regression_parameters_11776.20080105.txt';
rfilename = 'Residuals_11776.20080105.txt';

%%% Read data from input files

headerformat = repmat('%s ', 1, 2 + ns);
dataformat = [ repmat('%s ', 1, 2) repmat('%f ', 1, ns) ];

fid = fopen(cfilename);
cheader = textscan(fid, headerformat, 1, 'delimiter', '\t');
ctable = textscan(fid, dataformat, ng, 'delimiter', '\t');
fclose(fid);

fid = fopen(efilename);
eheader = textscan(fid, headerformat, 1, 'delimiter', '\t');
etable = textscan(fid, dataformat, ng, 'delimiter', '\t');
fclose(fid);

%%% Re-shape data and labels into more conveniently manipulated forms

c = zeros(ng,ns);
e = zeros(ng,ns);
sample_label{ns} = [];
gene_label{ng} = [];
gene_name{ng} = [];

for s = 1:ns
    c(:,s) = ctable{1,s+2}(:);
    e(:,s) = etable{1,s+2}(:);
    sample_label{s} = cheader{1,s+2}{1};
end

for g = 1:ng
    gene_label{g} = ctable{1,1}{g};
    gene_name{g} = ctable{1,2}{g};
end

clear ctable etable cheader eheader;

%%% Next: deal with NaNs
%
% In the given data set, there are 500 genes for which
% all 80 samples have copy-number value of NaN.
%
% This script handles the more general case where NaNs
% may occur more sporadically throughout the data.
% Instead of discarding any gene that has at least one NaN,
% NaNs are treated as missing data points, and the
% number (n) of data points in the regression calculation is
% decreased accordingly.  Genes with no good data points (n==0)
% will invoke DIV0 during regression and give a,b,r = NaN

%%% Count good data points (n) for each gene

good = ~isnan(c) & ~isnan(e);
n = sum(good,2);

%%% Transfer data to x and y, replacing bad data points with (0,0)

x = c;
y = e;

x(~good) = 0;
y(~good) = 0;

%%% Perform linear regression
%
%     given:
%         x = copy number
%         y = gene expression
%         n = number of good data points, for each gene
%
%     compute:
%         a = slope
%         b = intercept
%         r = Pearson correlation

xx = x.^2;
yy = y.^2;
xy = x.*y;

sx = sum(x,2);
sy = sum(y,2);
sxx = sum(xx,2);
syy = sum(yy,2);
sxy = sum(xy,2);

ntor = (n.*sxy)-(sx.*sy);
r = ntor ./ (sqrt((n.*sxx)-(sx.^2)) .* sqrt((n.*syy)-(sy.^2)));
a = ntor ./ ((n.*sxx)-(sx.^2));
b = (sy-(a.*sx)) ./ n;

%%% Alternatively, polyfit can be used to do the linear regresion.
%%% This gives exactly the same results for genes having no NaNs.
%%% For genes having a NaN anywhere in c or e, polyfit gives a,b = NaN
% 
%   a = zeros(ng,1);
%   b = zeros(ng,1);
% 
%   for g=1:ng
%       p = polyfit(c(g,:), e(g,:), 1);
%       a(g) = p(1);
%       b(g) = p(2);
%   end

%%% Compute residuals

efit = ((a*ones(1,ns)) .* c) + (b*ones(1,ns));
resid = e - efit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Write results to output files
%%%

%%% Output regression parameters

fid = fopen(pfilename, 'w');
fprintf(fid, 'REGRESSION\tNAME\ta\tb\tr\n');
for g=1:ng
    fprintf(fid, '%s\t%s\t%f\t%f\t%f\n', ...
        gene_label{g}, gene_name{g}, a(g), b(g), r(g));
end
fclose(fid);

%%% Output residuals

fid = fopen(rfilename, 'w');
headerformat = [repmat('%s\t', 1, 2+ns) '\n'];
fprintf(fid, headerformat, 'RESIDUALS', 'NAME', sample_label{:});

dataformat = ['%s\t%s\t' repmat('%f\t', 1, ns) '\n'];
for g=1:ng
    fprintf(fid, dataformat, gene_label{g}, gene_name{g}, resid(g,:));
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Generate figures
%%%

set(0, 'Units', 'pixels');
scnsize = get(0, 'ScreenSize');
w = scnsize(3);
h = scnsize(4);

% Fig. 1: Histogram of correlation coefficients

figure(1);
set(gcf, 'position', [0.1*w 0.25*h 0.75*h 0.6*h]);
hist(r,20);
title('Correlation of gene expression to copy number');
xlabel('Pearson correlation coefficient R');
ylabel('number of genes');
g = findobj(gca, 'Type','patch');
set(g,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k');
annotation(figure(1),'line',[0.385 0.385],[0.106 0.92],...
    'LineWidth',5,'Color',[1 0 0]);

% Fig. 2: Bitmap representation of c

% depict copy number as        c        cdraw
%                            -1.0 blue     1
%                             0.0 white  101
%                            +1.0 red    201
%                             NaN grey   202

cdraw = max(-1,min(c,1));
cdraw = (cdraw*100) + 101;
cdraw(isnan(c)) = 202;

cmap = ones(211,3);
cmap(1:101,1) = 0:0.01:1;
cmap(1:101,2) = 0:0.01:1;
cmap(101:201,2) = 1:-0.01:0;
cmap(101:201,3) = 1:-0.01:0;
cmap(202:211,:) = repmat([0.8 0.8 0.8],10,1);

figure(2);
set(gcf, 'position', [0.04*w 0.05*h, 0.28*w, 0.6*h]);
colormap(cmap);
image(cdraw);
title('copy number');
xlabel('samples');
ylabel('genes');
colorbar('YTick', [6 101 196 207], 'YTickLabel', ...
     {'Deleted','Neutral','Amplif.', 'No data'});

% Fig. 3: Bitmap representation of e

% depict expression data as     e        edraw
%                               4          1
%                              15         64
% using default MATLAB colormap

edraw = max(4,min(e,15));
edraw = ((edraw-4) * (63/11)) + 1;

figure(3)
set(gcf, 'position', [0.36*w 0.05*h, 0.28*w, 0.6*h]);
colormap('default');
image(edraw);
title('gene expression');
xlabel('samples');
ylabel('genes');
colorbar('Ytick', []);

% Fig. 4: Bitmap representation of resid

% depict residuals as       resid     residdraw
%                             -1         1
%                          0,NaN        32.5
%                             +1        64
% using default MATLAB colormap

residdraw = max(-1,min(resid,1));
residdraw(isnan(resid)) = 0;
residdraw = ((residdraw+1) * (63/2)) + 1;

figure(4);
set(gcf, 'position', [0.68*w 0.05*h, 0.28*w, 0.6*h]);
colormap('default');
image(residdraw);
title('residuals');
xlabel('samples');
ylabel('genes');
colorbar('Ytick', []);

% finally, bring Fig.1 to foreground

figure(1);

