function heatmap3d(D,cn_col,supid,smooth,pad,threshold)

figure(1); clf;
X=D.dat;

if exist('cn_col','var') && cn_col==1
  Y=dna_norm(X);
  Y(Y<-3)=-3;
  Y(Y>3)=3;
end


if exist('smooth','var')
  % smooth=1/sqrt(2);
  st=floor(smooth*3);
  sz=2*st+1;
  x=repmat([-st:1:st],sz,1);
  y=x';
  B=exp(-(((x/smooth).^2)/2+((y/smooth).^2)/2));
  B=B./sum(B(:));
else
  B=1;
  st=0;
end
if ~exist('pad','var')
  pad=0;
end
XX=pad*ones(size(X,1)+2*st,size(X,2)+2*st);
XX((st+1):(end-st),(st+1):(end-st))=X;
A=conv2(XX,B,'same');
A=A((st+1):(end-st),(st+1):(end-st));
if exist('threshold','var')
  A(A<threshold)=threshold;
end
%s=surf(1:size(A,1),1:size(A,2),A(:),Y(:));
if ~exist('Y','var')
  Y=A;
end
s=surf(1:size(A,2),1:size(A,1),A,Y);
set(s,'LineStyle','none')
shading interp
material shiny
% set(gca,'CameraPosition',[-235.397 -1014.41 99.4226] ,'CameraUpVector',[1.895 7.60998 0.743145]);
axis tight
ax=axis;
ysz=ax(4)-ax(3);
zsz=ax(6)-ax(5);
if exist('supid','var')
  for i=1:length(supid)
    supid_i=supid(i);
    ax=axis;
    axis([ax(1:2) ax(3)-ysz*(D.supmark(supid_i).y_3d+D.supmark(supid_i).y_3d_gap) ax(4) ax(5:6)]); 
    c=D.supmark(supid_i).colormap;
    v=D.supdat(supid_i,:);
    pw=D.supmark(supid_i).x_3d;
    r=0;
    for i=1:length(v)
      p{i}=draw_block([r+i-pw/2 r+i+pw/2],...
                      [ax(3)-ysz*(D.supmark(supid_i).y_3d+D.supmark(supid_i).y_3d_gap) ...
                       ax(3)-ysz*D.supmark(supid_i).y_3d_gap],[ax(5) ax(5)+D.supmark(supid_i).z_3d*zsz],c(v(i),:));
      set(p{i},'LineStyle','none');
    end
  end
end

%camproj('perspective')
axis off
bluepink

%m=mean(XX(~isnan(XX)));
%s=std(XX(~isnan(XX)));
%c=[m-3*s m+3*s];
%caxis(c);
