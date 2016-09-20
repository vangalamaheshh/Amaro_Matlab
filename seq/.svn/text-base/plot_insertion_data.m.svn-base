function plot_insertion_data(sample,chr,pos,win,t,n,sw,sw2)
if ~exist('sw','var'),sw = 300;end
if ~exist('sw2','var'),sw2 = 50;end

len = length(t.fhm_sm);
figure(1);clf;hold on;
set(gca,'ColorOrder',[0 0 0;1 0 0;0 0 0; 1 0 0;0 0 1;0 0 1;0 0.7 0.2; 0 0.7 0.2]);
sc = 50; o = -0.5;
plot([t.fhm_sm t.rhm_sm o+n.fhm_sm o+n.rhm_sm t.dhmz/sc o+n.dhmz/sc t.shmz/sc o+n.shmz/sc]);
xlim([sw len-sw]); ylim([o-0.3 o+1]);set(gca,'position',[0.02 0.05 0.96 0.9]);
text(0.05*len-2*sw,o+0.9,[sample '  chr' num2str(chr) ':' num2str(pos) ' +/- ' num2str(win)],'fontsize',15);
text(0.05*len-2*sw,o+0.8,['smooth\_window = ' num2str(sw)],'fontsize',12);
text(0.05*len-2*sw,0.15,'T');text(0.05*len-2*sw,o+0.15,'N');
hold off
