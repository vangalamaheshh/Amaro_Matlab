function [b,pks,bootsam]=bootstrap_peak_regions(C,nboot,score_type,qv_thresh,ext,nparts,lsfdir)

addpath ~gadgetz/matlab/snp

if ischar(C)
  Cname=C;
  tmp=load(Cname);
  C=tmp.C;
end

nsample=size(C.dat,2);

if exist('nparts','var') && nparts>1
  l=lsf(lsfdir);
  bidx=get_parts(1:nboot,nparts);
  h=zeros(nparts,1);
  if exist('Cname','var')
    C=Cname;
  end
  for i=1:nparts
    [l,h(i)]=bsub(l,{'b','pks','bootsam'},'bootstrap_peak_regions',...
                  {C,length(bidx{i}),score_type,qv_thresh,ext});
  end
  [l,res]=wait(l); % wait for all
  pks=[];
  bootsam=zeros(nsample,nboot);
  for k=1:2
    b{k}=zeros(size(res{h(1)}.b{k},1),1);
  end
  for i=1:nparts
    pks=[ pks res{h(i)}.pks];
    bootsam(:,bidx{i})=res{h(i)}.bootsam;
    for k=1:2
      b{k}=b{k}+res{h(i)}.b{k};
    end
  end    
else
  pks=[];
  bootsam=ceil(nsample*rand(nsample,nboot));
  for bi=1:nboot
    sam=bootsam(:,bi);
    CL=reorder_D_cols(C,unique(sam)); % unique or not unique?
    [regs,pvs]=bootstrap_peak_regions_step(CL,score_type,qv_thresh,ext);
    pks(bi).pvs=pvs;
    for k=1:2
      tmp=rmfield(regs{k},'segments');
      regs{k}=tmp;
    end
    pks(bi).regs=regs;
  end

  for k=1:2
    b{k}=zeros(size(C.dat,1),1);
    for i=1:length(pks)
      for j=1:length(pks(i).regs{k})
        b{k}(pks(i).regs{k}(j).peak_st:pks(i).regs{k}(j).peak_en)=...
            b{k}(pks(i).regs{k}(j).peak_st:pks(i).regs{k}(j).peak_en)+1;
      end
    end
  end
end




