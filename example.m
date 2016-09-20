%This script gives scatter plot, contour plot and ridgelines based on HMAC
%and REM. This plots are available only for 2-dim data. 

%import data
%Here 'glass2d.txt' is used for example.
data=textread('glass2d.txt');
s=max(std(data));sigmas=s*0.1:s*1.9/19:s*2;
[n_cluster,level, modes, members]=HMAC(data, sigmas);
i=2; %choose cluster level to get ridgelines and separability
k=min(find(level==i));
mode=eval(['modes.c',int2str(i)]);
[ridge,S,sep,density]=REM(data, mode,sigmas(k), members(:,i), 0:0.05:1);

%draw scatter plots
subplot(1,2,1);gscatter(data(:,1),data(:,2),members(:,i),'brgmck','.ox+*sdv^<>ph')
xaxis=xlim;yaxis=ylim; %save limits of axis for contour plot

%draw contour plots
[X,Y] = meshgrid(xaxis(1):.1:xaxis(2), yaxis(1):.1:yaxis(2)); 
Xn=reshape(X,prod(size(X)),1);Yn=reshape(Y,prod(size(Y)),1);
h=zeros(size(Xn));
for i=1:length(data)
    h=h+1/length(data)*mvnpdf([Xn,Yn],data(i,:),diag(sigmas(k)*ones(1,2)));
end
subplot(1,2,2);contour(X,Y,reshape(h,size(X)),100)

%add ridge lines to the contour plot
hold on
i=1; j=2; %choose two clusters which will be connected by a ridgeline: i<j 
ridge=eval(['ridge.c',int2str(i),'c',int2str(j)]);
plot(ridge(:,1),ridge(:,2))
plot(mode(:,1),mode(:,2),'x') %point out centers with symbol 'x'
hold off
          