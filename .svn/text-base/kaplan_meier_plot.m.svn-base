function [T,S]=kaplan_meier_plot(t,c,color,pch)
v=t;
vnotc=v(find(c==0));
t=unique(vnotc);
T=[];
S=[];
survival=1;
var_h=0;
ts=[0];
sur=[1];
b=1;
T(1)=0;
S(1)=1;

for i=1:length(t)
  ti=t(i);
  nrisk=length(find(v>=ti));
  nevent=length(find(vnotc==ti));
  survival=survival*(1-nevent/nrisk); % Kaplan-Meier
  var_h=var_h+nevent/(nrisk*(nrisk-nevent)); % Greenwood
					     % var_h=var_h+nevent*(nrisk-nevent)/(nrisk^3); % Klein 
  stderr=sqrt(survival^2*var_h);
  z=norminv(1-0.05/2);
  h=-log(survival);
  ub=exp(-h*exp(-z*sqrt(var_h)/h));
  lb=exp(-h*exp(+z*sqrt(var_h)/h));
  
%  ub=exp(-exp(log(h+z*sqrt(var_h)/h )));
%  lb=exp(-exp(log(h-z*sqrt(var_h)/h )));
%     disp([ti nrisk nevent survival stderr h sqrt(var_h) ...
%             h+z*sqrt(var_h)/h h-z*sqrt(var_h)/h]);
% disp([ti survival ]);
  b=b+2;
  T(b-1)=ti;
  S(b-1)=S(b-2);
  T(b)=ti;
  S(b)=survival;
  UB(b)=ub;
  LB(b)=lb;
  UB(b-1)=UB(b-2);
  LB(b-1)=LB(b-2);
end
%keyboard
S(b+1)=S(b);
UB(b+1)=UB(b);
LB(b+1)=LB(b);
T(b+1)=max(v)+1;
H=figure(1);
hold on;
if pch
    for i=1:2:length(T)
        if (S(i)<1)
            ph=patch([ T(i) T(i+1) T(i+1) T(i) T(i)], ...
                [LB(i) LB(i) UB(i) UB(i) LB(i)],...
                [0.7 0.7 0.7]);
            set(ph,'LineStyle','none','FaceAlpha',0.5);
        end
    end
end
plot(T,S,color,'LineWidth',1);

%errorbar(T,S,S-LB,UB-S);

ind_c=find(c);
for i=1:length(ind_c)
  T_c=v(ind_c(i));
  m=max(T(find(T<v(ind_c(i)))));
  S_c=min(S(find(ismember(T,m))));
  plot(T_c,S_c,[ color '+'],'MarkerSize',10,'LineWidth',1);
end

a=axis;
axis([a(1) a(2) 0 a(4)]);




