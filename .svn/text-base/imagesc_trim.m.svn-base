function imagesc_trim(D,c,robust_flag)
% imagesc_trim(D)
%    performs IMAGESC and sets the color scale to span 3 standard 
%    deivations below and above the mean. 
%
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%


if exist('robust_flag','var') && robust_flag
  m=median(D(~isnan(D)));
  s=mad(D(~isnan(D)));
else
  m=mean(D(~isnan(D)));
  s=std(D(~isnan(D)));
end

Dn=isnan(D);
nanpos=find(Dn);

if ~exist('c','var') || isempty(c)
    c=[m-3*s m+3*s];
end

if ~isempty(nanpos)
    Dx=floor((D-c(1))/(c(2)-c(1)+eps)*64)+1;
    Dx(Dx<1)=1;
    Dx(Dx>64)=64;
    Dx(nanpos)=1;
    cm=colormap;
    Dc=cm(Dx,:);
    Dc(nanpos,:)=repmat([90 100 90]/255,length(nanpos),1);
    Dc=reshape(Dc,size(Dx,1),size(Dx,2),3);
    image(Dc);
else
    imagesc(D);
    caxis(c);
end
