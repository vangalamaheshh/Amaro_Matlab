function jitterplot(datamatrix,colheaders,uselog,order,outdir,jittervalue,clrs,dotsize,drawmedians,medianclr)
% draw vertical strip-plots with jitter (like boxplot outliers, without the
%   boxplot). Only the first input is strictly required.
%
% datamatrix = 2-D matrix of data, MxN where M is data for each sample, and
%   N is the number of columns
% colheaders = cell array of N column headers (generate with import()).
% uselog = 0 for linear, 1 for log Y-scale
% order = vector of N integers indicating a new location on which to
%   display the column. (default is [1:length(colheaders)]). For example,
%   if columns are [A B C D E], order [1 3 4 5 2] would display in the new 
%   order [A E B C D].
% outdir = directory to which to write jitterplot.pdf (not currently
%   enabled)
% jittervalue = number between 0 and 1, determines the width of the jitter
% clrs = 1x3 scalar color value, (default [0 0 0]=black)
% dotsize = scalar pixels^2 area of dots size (default 25)
% drawmedians = 0 for no, 1 for yes (default 1)
% medianclr = color for median line (default [0 0 1]=blue)
% 
% Peleg Horowitz, Broad Institute, 2012

% close all
dbstop if error

if ~exist('datamatrix','var')
    error('must provide data!')
end

colnum=size(datamatrix,2);
rownum=size(datamatrix,1);

if ~exist('colheaders','var') || isempty(colheaders)
    colheaders=cell(1,colnum);
    for i=1:colnum
        colheaders{i}='';
    end
end

if size(datamatrix,2) ~= size(colheaders,2)
    error('number of columns and column headers do not match!');
end

if ~exist('order','var') || isempty(order)
    order=[1:colnum];
else
    if length(order)~=colnum
        error('length of order field does not match number of columns!');
    end
end

if ~exist('outdir','var') || isempty(outdir)
    outdir = pwd;              
end
%add_slash_if_needed(outdir);

if ~exist('jittervalue','var') || isempty(jittervalue)
    jittervalue=0.5;
end

if ~exist('uselog','var') || isempty(uselog)
    uselog=0;
end

if ~exist('clrs','var') || isempty(clrs)
    clrs=[0 0 0];
end

if ~exist('dotsize','var') || isempty(dotsize)
    dotsize=25;
end

if ~exist('drawmedians','var') || isempty(drawmedians)
    drawmedians=1;
end

if ~exist('medianclr','var') || isempty(medianclr)
    medianclr=[0 0 1];
end

xs=nan(colnum*rownum,1);
ys=xs;
med=nan(colnum,1);

for i=1:colnum
    ii=order(i);
    ys((i-1)*rownum+1:i*rownum)=datamatrix(:,i);
    xs((i-1)*rownum+1:i*rownum)=jittervalue*rand(rownum,1)+ii-0.5-jittervalue/2;
    med(ii)=nanmedian(datamatrix(:,i));
end

xs=xs(~isnan(ys));
ys=ys(~isnan(ys));

if uselog % fix -infinity for log(0) scores
    if min(ys)==0
        floor=min(ys(ys~=0));
        for i=1:length(ys)
            if ys(i)==0
                ys(i)=ys(i)+floor/10;
            end
        end
    end
end

% scatter plot
figure
scatter(xs,ys,dotsize,clrs);
hold on

if drawmedians
    for i=1:colnum
        wd=max(jittervalue+0.2,0.3);
        line(i-0.5+[-wd/2 wd/2],[med(i) med(i)],'Color',medianclr);
    end
end

tmp=colheaders;
colheaders(order)=tmp(:);

set(gca,'XTick',0.5:1:colnum-0.5);
if exist('colheaders','var') && ~isempty(colheaders)
    set(gca,'XTickLabel',colheaders);
end

if uselog
    set(gca,'YScale','log')
end




