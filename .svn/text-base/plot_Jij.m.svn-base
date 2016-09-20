function plot_Jij(dat,Jij)

clf;
set(gcf,'renderer','zbuffer');
hp=plot(dat(:,1),dat(:,2),'x');

for i=1:size(Jij,1)
  hl(i)=line([dat(Jij(i,1),1) dat(Jij(i,2),1)],[dat(Jij(i,1),2) ...
                      dat(Jij(i,2),2)],'Color','red');
  set(hl(i),'Tag',[num2str(i) ' ' num2str(Jij(i,1)) ':' num2str(Jij(i,2))]);
end
