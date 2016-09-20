prop_p=0.1;
prop_c=0.3;
prop_pc=0.1; % P(p|c)
M=4000;
M_pc=[0 M*prop_c*(1-prop_pc) M*prop_p-M*prop_c*prop_pc M*prop_c* ...
      prop_pc];
M_pc(1)=M-sum(M_pc);

delta_p=1;
delta_pc=0;
sig=1;
delta_c=1;

D=sim_data(M_pc,20,0.7,sig,delta_p,delta_c,delta_pc);
D=reorder_D_cols(D,1:(size(D.dat,2)-3));


