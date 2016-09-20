function plot_broad_results(out_dir,D,rg_file,qvthresh,ignore_XY,sigcap,sample)
% A function to run at the end of GISTIC 2.0 in order to graphically plot
% the broad significance results. Requires that GISTIC write a file, 
% broad_results.mat, into the output directory (out_dir).
% The q-values are plotted on log scale (unlike skyline focal gistic plots 
% which are on log-log scale). For plotting purposes, q-values are capped 
% at 10^-8.
%
% Parameters:
% 
%     out_dir: the GISTIC output directory in which the broad_results.mat
%              file has been written (string)
%     D:       the D struct at the end of GISTIC; only the D.chrn field is 
%              used to determine max chromosome number and lengths (struct)
%              - can pass through a field D.plotY set to 0 to suppress
%              plotting of Y chromosome.
%     rg_file: the name/location of the .mat reference genome file (string)
%              from which relative arm lengths are determined for p and q
%    qvthresh: qv threshold, defaults to 0.25
%   ignore_xy: defaults to 0; otherwise, excludes X and Y
%      sigcap: caps plot to maximum of sigcap, defaults to 10^-30
%      sample: sampling rate for large figure (to reduce resolution)
%
%
% Peleg Horowitz, 2012

if ~exist('D.plotY','var')
    D.plotY=0;
end

if ~exist('qvthresh','var')
    qvthresh=0.25;
end
log_qvthresh = -log10(qvthresh);

if ~exist('ignore_XY','var')
    ignore_XY = 0;
end

if~exist('sigcap','var')
    sigcap = 30;
end

if ~exist('sample','var')
    sample = 1;
end

% lay out figure sub-panel dimensions for later use
lt = [0.1 0.1 0.37 0.8];
ct = [0.475 0.1 0.05 0.8];
rt = [0.53 0.1 0.37 0.8];

% load broad data

out_dir=add_slash_if_needed(out_dir);
brod = load([out_dir, 'broad_results.mat']);

D = reorder_D_rows(D,[1:sample:size(D.dat,1)]);
brod = brod(1:sample:end);

load(rg_file);
brod.logqA = -log10(brod.qA);
brod.logqD = -log10(brod.qD);

% find the integer above the highest -log10(q) value
% cap the log score at 10^8 if needed
maxbrod = round(max(max(brod.logqA),max(brod.logqD))+0.5);
maxdel = max(brod.logqD);
maxamp = max(brod.logqA);

if maxbrod > sigcap
    maxbrod = sigcap;
    brod.logqA = min(brod.logqA,sigcap);
    brod.logqD = min(brod.logqD,sigcap);
    if maxamp>sigcap
        capped{1}='<';
    else
        capped{1}='';
    end
    if maxdel>sigcap
        capped{2}='<';
    else
        capped{2}='';
    end
else
    capped{1}='';
    capped{2}='';
end

if ignore_XY==0
    maxchrom = max(D.chrn);
else
    maxchrom = 22;
end

if maxchrom==24 && ~D.plotY
    maxchrom=23;
end
brod.chrn = zeros(maxchrom*2,1); brod.armn = zeros(maxchrom*2,1); brod.namen = cell(maxchrom*2,1);
brod.nlogqA = NaN(maxchrom*2,1); brod.nlogqD = NaN(maxchrom*2,1); brod.armlen = NaN(maxchrom*2,1);

for i=1:maxchrom % fill in zeros for the arms for which there is no data
    for j = 1:2
        brod.chrn(i*2-2+j)=i;
        brod.armn(i*2-2+j)=j-1;
        
        if j==1
            arm='p';
        else
            arm='q';
        end
        
        if i<23
            brod.namen{i*2-2+j}=[num2str(i) arm];
        elseif i==23
            brod.namen{i*2-2+j}=['X' arm];
        elseif i==24
            brod.namen{i*2-2+j}=['Y' arm];
        end
        
        if sum(strcmp([num2str(i) arm],brod.names))
            brod.nlogqA(i*2-2+j)=brod.logqA(strcmp([num2str(i) arm],brod.names));
            brod.nlogqD(i*2-2+j)=brod.logqD(strcmp([num2str(i) arm],brod.names));
        else
            brod.nlogqA(i*2-2+j)=0;
            brod.nlogqD(i*2-2+j)=0;
        end
    end
end

for i=1:maxchrom % set the arm lengths for chromosomes
    ps=length(grep(brod.namen(i*2-1),{cyto.name},1));
    qs=length(grep(brod.namen(i*2),{cyto.name},1));
    psandqs=ps+qs;
    chrlength=sum(D.chrn==i);
    brod.armlen(i*2-1) = ps/psandqs*chrlength;
    brod.armlen(i*2) = qs/psandqs*chrlength;
end

genomelen=sum(brod.armlen);
brod.armlen=brod.armlen/genomelen*maxbrod*lt(4)/lt(3); % fix the axes ratio to 1:1 by correcting for genome length then aspect ratio of the panel

% initialize nan matrices of x, y coordinates for the bars corners
x1=NaN(maxchrom*2,1);
x2a=x1; x2d=x1; xmax=x1; y1=x1; y2=x1;
ampx=NaN(4,maxchrom*2);
delx=ampx; bkgx=ampx; ydat=ampx;

% fill in the matrix with q and armlength dimensions for the bars
for i=1:maxchrom*2
    x1(i)=0;
    x2a(i)=brod.nlogqA(i);
    x2d(i)=brod.nlogqD(i);
    xmax(i)=maxbrod;
    y1(i)=sum(brod.armlen(1:end))-(sum(brod.armlen(1:i))-brod.armlen(i));
    y2(i)=sum(brod.armlen(1:end))-sum(brod.armlen(1:i));
    ampx(:,i)=[x1(i); x2a(i); x2a(i); x1(i)];
    delx(:,i)=[maxbrod-x1(i); maxbrod-x2d(i); maxbrod-x2d(i); maxbrod-x1(i)];
    bkgx(:,i)=[x1(i); xmax(i); xmax(i); x1(i)];
    ydat(:,i)=[y1(i); y1(i); y2(i); y2(i)];
end

% initialize figure and subpanels
figure; clf;
set(gcf,'Color',[1 1 1]);
set(gcf,'Name','Broad GISTIC analysis');
hlt = subplot('Position',lt);
hct = subplot('Position',ct);
hrt = subplot('Position',rt);
set(hct,'Xtick',[],'Ytick',[],'Box','off','XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
set(hlt,'XLim',[0 maxbrod],'YLim',[0 sum(brod.armlen(1:end))]);
set(hrt,'XLim',[0 maxbrod],'YLim',[0 sum(brod.armlen(1:end))]);

darkgreen=[91 186 71]/255;

% plot the bars for amps (1) and dels (2)
for k=1:2
    if k==1
        xdat=ampx;
        color = [1 0 0]; % red for amps
        subplot(hrt);
        titl = 'Gains';
    else
        xdat=delx;
        color = [0 0 1]; % blue for dels
        subplot(hlt);
        titl = 'Losses';
    end
    
    % set up tick marks and tick mark labels for the X-axis
    if maxbrod >15
        tickmarks=[0 round(maxbrod/4) round(maxbrod/2) round(3*maxbrod/4)  maxbrod];
        if k==2
            tickmarks=sort(maxbrod-tickmarks);
            ticklabels={[capped{k} '10^-' num2str(maxbrod)] ['10^-' num2str(tickmarks(4))] ['10^-' num2str(tickmarks(3))] ['10^-' num2str(tickmarks(2))] '1'};
        else
            ticklabels={'1' ['10^-' num2str(tickmarks(2))] ['10^-'  num2str(tickmarks(3))] ['10^-' num2str(tickmarks(4))] [capped{k} '10^-' num2str(maxbrod)]};
        end        
    elseif maxbrod >8
        tickmarks=[0 1 3 6 maxbrod];
        
        if k==2
            tickmarks=sort(maxbrod-tickmarks);
            ticklabels={[capped{k} '10^-' num2str(maxbrod)] '10^-6' '10^-3' '0.1' '1'};
        else
            ticklabels={'1' '0.1' '10^-3' '10^-6' [capped{k} '10^-' num2str(maxbrod)]};
        end        
    elseif maxbrod>5 %6 to 8
        tickmarks=[0 1 2 4 maxbrod];
        
        if k==2
            tickmarks=sort(maxbrod-tickmarks);
            ticklabels={[capped{k} '10^-' num2str(maxbrod)] '10^-4' '0.01' '0.1' '1'};
        else
            ticklabels={'1' '0.1' '0.01' '10^-4' [capped{k} '10^-' num2str(maxbrod)]};
        end
    elseif maxbrod>2 %3 to 5
        tickmarks=[0 1 2 maxbrod];
        
        if k==2
            tickmarks=sort(maxbrod-tickmarks);
            ticklabels=[10^(-maxbrod) 0.01 0.1 1];
        else
            ticklabels=[1 0.1 0.01 10^(-maxbrod)];
        end
    else % <3
        tickmarks=[0 1 2];
        
        if k==2
            tickmarks=sort(maxbrod-tickmarks);
            ticklabels=[0.01 0.1 1];
        else
            ticklabels=[1 0.1 0.01];
        end
    end
         
    set(gca, 'XTick', tickmarks, 'XTickLabel', ticklabels);
    set(gca, 'Ytick', [], 'Ztick', []);
    set(get(gca, 'Title'), 'String', titl);
    set(get(gca, 'XLabel'), 'String', 'FDR q-value');
            
    bkgcol = 1-mod(brod.chrn+1,2)*0.1; %alternating grey-white background for chromosomes
    
    for i=1:maxchrom*2
        phbkg = patch(bkgx(:,i),ydat(:,i),[bkgcol(i,:) bkgcol(i,:) bkgcol(i,:)]);
        set(phbkg,'FaceAlpha',0.9,'EdgeColor','none');
        if mod(i,2) % dashed line at centromeres
            lh=line([0 maxbrod],[sum(brod.armlen(1:end))-sum(brod.armlen(1:i)) sum(brod.armlen(1:end))-sum(brod.armlen(1:i))],'LineStyle',':');
            set(lh,'Color',[0.3 0.3 0.3]);
        end
    end
        
    ph = patch(xdat,ydat,color);
    set(ph,'FaceAlpha',0.9,'EdgeColor','none');
    
    if k==1 % draw q-value significance threshold 
        line([log_qvthresh log_qvthresh],[0 sum(brod.armlen(1:end))],'Color',darkgreen);
    else
        line([maxbrod-log_qvthresh maxbrod-log_qvthresh],[0 sum(brod.armlen(1:end))],'Color',darkgreen);
    end
    
    box on
end

% plot the genome legend
set(hct,'XLim',[0 100]);
set(hct,'YLim',[0 sum(brod.armlen(1:end))]);
subplot(hct);

x1=NaN(maxchrom,1); % initialize empty x-y matrices
x2=x1; xt=x1; y1=x1; y2=x1; y3=x1; y4=x1; yt=x1; 
chrx=NaN(4,maxchrom); chry=chrx; ceny=chrx;

for i=1:maxchrom % fill in values for x,y coordinates of each corner of each chromosome
    x1(i)=22+mod(i,2)*40;
    x2(i)=x1(i)+16;
    xt(i)=30+mod(i+1,2)*40;
    y1(i)=sum(brod.armlen(1:end))-(sum(brod.armlen(1:i*2))-brod.armlen(i*2)-brod.armlen(i*2-1));
    y2(i)=sum(brod.armlen(1:end))-sum(brod.armlen(1:i*2));
    y3(i)=sum(brod.armlen(1:end))-sum(brod.armlen(1:i*2-1));
    y4(i)=min(y3(i)+sum(brod.armlen(1:end))*0.005,sum(brod.armlen(1:end)));
    yt(i)=(y1(i)+y2(i))/2;
    chrx(:,i)=[x1(i) x2(i) x2(i) x1(i)];
    chry(:,i)=[y1(i) y1(i) y2(i) y2(i)];
    ceny(:,i)=[y3(i) y3(i) y4(i) y4(i)];
end

phc=patch(chrx,chry,[0.8 0.8 0.8]); 
phcc=patch(chrx,ceny,[0.6 0.6 0.6]);
set(phc,'FaceAlpha',0.9,'EdgeColor','none');
set(phcc,'FaceAlpha',0.9,'EdgeColor','none');

%{
for i=1:maxchrom    %add numbers of chromosomes
    zt=num2str(i);
    if i==23
        zt='X';
    elseif i==24
        zt='';  % Y is too small to write next to
    end
    
    handl=text(xt(i),yt(i),zt,'HorizontalAlignment','center');
    if i<20
        set(handl,'FontSize',10);
    elseif i<23
        set(handl,'FontSize',9);
    else % X
        set(handl,'FontSize',10);
    end
end
%}

fprintf('Done!\n');

% write figure file
pathname = [out_dir 'broad_significance'];
saveas(gcf,[pathname '.fig'],'fig');
print -painters -dpdf [pathname '.pdf'];
saveas(gcf,[pathname '.png'],'png');

