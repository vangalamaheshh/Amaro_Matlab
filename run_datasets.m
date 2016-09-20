cd ~/projects/booster/datasets
files={'duke_breast.mat','all_mll.mat','rosetta_breast_out.10K.mat','dlbcl_new_out.mat',...
        'dlbcl_media_new.mat','sim.mat'};
minvar=(0.2)^2;
gamma=95;

l=lsf('/xchip/data/gadgetz/lsfres/');
h=zeros(length(files),1);
for i=1:length(files)
  [l,h(i)]=bsub(l,{},'run_dataset',{files{i},minvar,gamma,0});
end
[l,res]=wait(l); % wait for all



% bsub -q normal -r -o run_stdout.txt -e run_stderr.txt matlab -nodisplay -r run_datasets; exit;
