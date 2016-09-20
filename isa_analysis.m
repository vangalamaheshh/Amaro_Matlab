function sets=isa_analysis(dat,sets)

% the data is already scaled and positive

dat=log2(dat);
dat=dat-repmat(mean(dat,2),1,size(dat,2));
dat(dat>5)=5;
dat(dat<-5)=-5;


% start around each gene


disp('start');
down=1;

EG=dna_norm(dat)./sqrt(size(dat,2)-1);
ECT=dna_norm(dat')'./sqrt(size(dat,1)-1);

EG=EG';

t_g=4; t_c=2;

disp('loop');
for k=1:5:23
  t_g=1.7+k*0.1;
  
  for j=1:length(sets)
    if mod(j,10)==0
      disp(j);
    end
    sg={}; sc={};
    tmp=sets{j}.u95';
    sg_tmp=sparse(size(EG,2),1);
    %    sg_tmp=zeros(size(EG,1),1);
    sg_tmp(tmp(find(~isnan(tmp))))=1;
    sg{1}=sg_tmp;
%    keyboard
    sc{1}=[];
    sets{j}.step=10;
    for i=2:10  
      [sg_,sc_]=isa_step(sg{i-1},t_g,t_c,EG,ECT,down); 
%      sg{i}=sparse(sg_);
%      sc{i}=sparse(sc_);
      sg{i}=sg_;
      sc{i}=sc_;
      old_sg=full(sg{i-1});
      if ~any( abs(sg_-old_sg) > max(repmat(0.01,length(sg_),1), ...
                                     0.01*abs(old_sg)) )
        sets{j}.step=i;
        break
      end
    end
    sets{j}.sg{k}=sg;
    sets{j}.sc{k}=sc;
  end
%  save results2.mat sets t_g t_c down k
end

