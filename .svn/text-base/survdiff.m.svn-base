function p=survdiff(v1,c1,v2,c2)

v1notc=v1(find(c1==0));
v2notc=v2(find(c2==0));
v12=sort([v1 v2]);
v12notc=sort([v1notc v2notc]);
t=unique([v1notc v2notc]); % T values for steps in the Kaplan-Meier curve
cup=0; cne1=0; cexp_ne1=0; cvar=0;
for i=1:length(t)
    ti=t(i);
    nrisk1=length(find(v1>=ti)); % number of patients from v1 with survival >= ti (which are at risk)
    if isempty(v1notc)
        nevent1=[];
    else
        nevent1=length(find(v1notc==ti)); % number of deaths from v1 (not including the censored) at ti 
    end
    nrisk12=length(find(v12>=ti)); % number of patients from v1 and v2 with survival >= ti (which are at risk)
    if isempty(v12notc)
        nevent12=[];
    else
        nevent12=length(find(v12notc==ti)); % number of deaths from v1 and v2 (not including the censored) at ti 
    end
    exp_ne1=nrisk1*(nevent12/nrisk12); % expected number of events in v1 according to the ratio in (v1 and v2)
    if (nrisk1==0) | (nrisk12==nrisk1) 
        var_ne1=0; 
    else
        var_ne1=exp_ne1*(nrisk12-nrisk1)*(nrisk12-nevent12)/(nrisk12*(nrisk12-1));
    end
    cne1=cne1+nevent1;
    cexp_ne1=cexp_ne1+exp_ne1;
    cvar=cvar+var_ne1;
    verbose([ ti nrisk1 nevent1 nrisk12 nevent12 exp_ne1 sqrt(var_ne1) cne1 cexp_ne1 (nevent1-exp_ne1)^2/var_ne1]);
%    disp([ ti nevent1 exp_ne1 sqrt(var_ne1) cne1 cexp_ne1 (nevent1-exp_ne1)^2/var_ne1 (cne1-cexp_ne1)^2/cvar]);
%    disp((cne1-cexp_ne1)^2/cvar);
end
p=1-chi2cdf((cne1-cexp_ne1)^2/cvar,1);


