function [bel,cor]=semisup_spc(J,params)

TrainLabels=params.TrainLabels;
n_train=size(TrainLabels,1);
if ~isfield(params,'q')
    q=max(TrainLabels(:,2));
else
    q=params.q;
end
% push points q forward and add interaction with labelled points
J(:,1:2)=J(:,1:2)+q;
J=sortrows([ TrainLabels(:,2) TrainLabels(:,1)+q repmat(1000,n_train,1); J]);

start_Js=min(find(J(:,1)>q));
Js=find(J(:,3)~=1000);
Jav=mean(J(Js,3));
N=max(max(J(:,[1 2])));

if ~isfield(params,'Ts')
    Ts=0.001:0.01:1.491
else
    Ts=params.Ts;
end 

p.n_bins=1000;
p.nrep=5;
p.cyc=300000;
p.final_cyc=50000;
p.nsample=1;
p.max_iter=1; %5
p.ignore_cyc_fraction=3000/p.cyc; %0.2;
p.slope=0;
p.wl_hist_ratio=0.8;
p.wl_min_hist_support=round(p.n_bins)*0.5;
p.final_f=0.01; % NOTE !!!!!
p.wl_f_start=0.01;
p.display=0;
p.wl_check_hist=2000;
p.wl_write_D=20000;

params=add_struct(p,params);


disp(sum(J(min(find(J(:,1)>q)):end,3)));
params.maxE=sum(J(min(find(J(:,1)>q)):end,3))*1;
%params.maxE=sum(J(min(find(J(:,1)>q)):end,3))*0.5; % 
params.beta=0;
params.alpha=1;

params.init_met=0;
params.init_maxE=params.maxE; 
params.init_beta=1/Ts(1);

M=length(Js);
Estep=(params.maxE-0)/(params.n_bins-1);
m=0:Estep:params.maxE;
m=m/Jav;

if isfield(params,'init_d')
    load(params.init_d)
else
    init_d=zeros(1,params.n_bins);
end 

has_labels=1;
obs_list={};
rand_seed=round(sum(100*clock));
disp(rand_seed);
rand('state',rand_seed);

params
D=[]; Obs={}; states=[]; states_obs=[];
belief=zeros(length(Ts),q*(N-q));
pair=zeros(length(Ts),size(J,1)-start_Js+1);
pairmat=zeros(length(Ts),(size(J,1)-start_Js+1)*q^2);
so={[],[],[],[]};
e_all=[];
Z=zeros(length(Ts),1);

lcf=zeros(length(Ts),1);
save_lcf=[];
save([ params.run_prefix '_BEFORE_RUN.mat']);
for k=1:params.nrep
  [D_k,Obs_k,states_k,states_obs_k]=potts_monte_carlo('WL',params,J,q,Ts,has_labels,obs_list,init_d);
  disp(['Finished rep ' num2str(k) ]);
  D=D_k;
  Obs={}; % not used;
%  states=[states states_k];
  states_obs=[states_obs states_obs_k];
   
  e=states_obs_k(1,:);
  params.minE=0;
  ebin=min(floor((e-params.minE)/Estep+0.5)+1,params.n_bins);
  e_hist=hist(ebin,1:params.n_bins);
  empty_e=find(e_hist==0);
  D(:,empty_e)=-Inf;
  d=D(end-1,:);
  n=N-q;
	
  state_prob=zeros(length(Ts),length(e));
  lcf_k=zeros(length(Ts),1);

%  ebin_hi=min(ceil((e-params.minE)/Estep)+1,params.n_bins);
%  ebin_low=ebin_hi-1;
%  ebin_frac=((e-params.minE)/Estep)-(ebin_low-1);
%  d_ebin_low=d(ebin_low);
%  d_ebin_hi=d(ebin_hi);
%  tmp=isinf(d_ebin_low);d_ebin_low(tmp)=d_ebin_hi(tmp);
%  tmp=isinf(d_ebin_hi);d_ebin_hi(tmp)=d_ebin_low(tmp);    
%  state_d=d_ebin_low.*(1-ebin_frac)+d_ebin_hi.*ebin_frac;

  state_d=d(ebin);

  for ti=1:length(Ts)
    T=Ts(ti);
    lcf_k(ti)=max(-e./T+state_d);
    p_e=exp(-e./T+state_d-lcf_k(ti));
    state_prob_k(ti,:)=p_e;
  end
  Z_k=sum(state_prob_k,2);
  
  save([params.run_prefix '_rep' sprintf('%d',k) '_before.mat'], ...
       'lcf_k','states_k','state_prob_k','Z_k','e','D','rand_seed');
  

  chunk_size=1000;
  i=1;
  belief_k=zeros(length(Ts),q*n);
  pair_k=zeros(length(Ts),size(J,1)-start_Js+1);
  pairmat_k=zeros(length(Ts),(size(J,1)-start_Js+1)*q^2);
  
  while (i <= length(e))
    range=i:(min(i+chunk_size-1,length(e)));
    chunk_state_belief=zeros(q*n,chunk_size);
    for j=1:q
      chunk_state_belief(((j-1)*n+1):(j*n),1:length(range))=states_k((q+1):end,range)==j;
    end
    
    chunk_state_pair=zeros(size(J,1)-start_Js+1,chunk_size);
    for j=1:size(chunk_state_pair,1)
      p_ij=J(j+start_Js-1,1:2);
      chunk_state_pair(j,1:length(range))=states_k(p_ij(1),range)==states_k(p_ij(2),range);
    end
    
    chunk_state_pairmat=zeros((size(J,1)-start_Js+1)*q^2,chunk_size);
    for j=1:((size(J,1)-start_Js+1))
      p_ij=J(j+start_Js-1,1:2);
      for q1=1:q
	for q2=1:q
	  chunk_state_pairmat((j-1)*q^2+(q1-1)*q+q2,1:length(range))=...
	      ((states_k(p_ij(1),range)==q1) & (states_k(p_ij(2),range)==q2));
	end
      end
    end
    
    belief_k=belief_k+state_prob_k(:,range)*chunk_state_belief(:,1:length(range))';
    pair_k=pair_k+state_prob_k(:,range)*chunk_state_pair(:,1:length(range))';
    pairmat_k=pairmat_k+state_prob_k(:,range)*chunk_state_pairmat(:,1:length(range))';
    i=i+length(range);
  end
  
  new_lcf=max([lcf lcf_k],[],2);
  f_old=exp(lcf-new_lcf);
  f_new=exp(lcf_k-new_lcf);
  save_lcf=[ save_lcf lcf_k ];

  belief=belief.*repmat(exp(f_old),1,size(belief,2))+belief_k.*repmat(exp(f_new),1,size(belief_k,2));
  pair=pair.*repmat(exp(f_old),1,size(pair,2))+pair_k.*repmat(exp(f_new),1,size(pair_k,2));
  pairmat=pairmat.*repmat(exp(f_old),1,size(pairmat,2))+pairmat_k.*repmat(exp(f_new),1,size(pairmat_k,2));
  Z=Z.*repmat(exp(f_old),1,size(k,2))+Z_k.*repmat(exp(f_new),1,size(Z_k,2));

  e_all=[e_all e];
  
  for i=1:4
    so_k{i}=sum(states_k==i);
    so{i}=[ so{i} so_k{i} ];
  end
  
  [se,si]=sort(e);
  states_k=states_k(:,si(1:10:end));

  save([params.run_prefix '_rep' sprintf('%d',k) '.mat'], ...
       'states_k','belief','pair','pairmat','Z','save_lcf','Z_k','belief_k','pair_k','pairmat_k','so_k','e','D','so','rand_seed');

  disp(['Finished calc ' num2str(k) ]);
end

bel=belief;
cor=pair;

