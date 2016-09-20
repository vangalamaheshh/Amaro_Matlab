function [e_amp,e_del]=generate_e(CL21);

xv=-4:0.01:4;
h=hist(CL21.dat(:),xv);
% plot(cumsum(h))
ph=(h+1);
ph=ph/sum(ph);
ph_amp=ph(xv>=0);
ph_amp(1)=sum(ph(xv<=0));
ph_amp=ph_amp/sum(ph_amp);
ph_del=ph(xv<=0);
ph_del(end)=sum(ph(xv>=0));
ph_del=ph_del/sum(ph_del);

e_amp=-log(ph_amp);
e_del=-log(ph_del);
e_del=fliplr(e_del);

e_amp=e_amp-e_amp(1);
e_del=e_del-e_del(1);

e_amp=e_amp/4;
e_del=e_del/4;
