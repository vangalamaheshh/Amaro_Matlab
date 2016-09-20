function [p] = script_for_A_p_sim_func(nSim)
tic
A=round(rand(nSim,1)*100);
p=sum(mod(A,10)~=2)/length(A)
toc
end