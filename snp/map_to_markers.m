function [C,q,ads,regs,map,map_back]=map_to_markers(CL21,q,ads,regs,Cnew,posfield)

if ~exist('posfield','var')
  pos=Cnew.pos;
else
  pos=getfield(Cnew,posfield);
end

CL21_save=CL21;
q_save=q;
ads_save=ads;
regs_save=regs;


emp=[];
map=zeros(size(Cnew.dat,1),1);
map_back=zeros(size(CL21.dat,1),1);
for i=1:23
  in_chr=find(CL21.chrn==i);
  in_chr_A=find(Cnew.chrn==i);
  for j=1:length(in_chr)
    if j==1
      st=0;
    else
      st=mean(CL21.pos(in_chr((j-1):j)));
    end
    
    if j==length(in_chr)
      en=Inf;
    else
      en=mean(CL21.pos(in_chr(j:(j+1))));
    end
    
    si=find(Cnew.chrn==i & pos>st & pos<en);
    if ~isempty(si)
      map(si)=in_chr(j);
    else
      [tmp,ti]=min(abs(pos(in_chr_A)-CL21.pos(in_chr(j))));
      map_back(in_chr(j))=ti+min(in_chr_A)-1;
      emp=[emp in_chr(j)];
    end
    if mod(j,100)==0
      disp(j);
    end
  end
  i
end

% save map_to_100K_cyto.mat map map_back

C.marker=Cnew.marker;
C.dat=zeros(size(Cnew.dat,1),size(CL21.dat,2));
C.pos=Cnew.pos;
C.chr=Cnew.chr;
C.chrn=Cnew.chrn;
C.dat=CL21.dat(map,:);

for k=1:2
  q{k}=ones(size(Cnew.dat,1),1);
  q{k}=q_save{k}(map);
  ads{k}=ones(size(Cnew.dat,1),1);
  ads{k}=ads_save{k}(map);
end
q1=q;
ads1=ads;
mrl=runlength(map_back);
for i=find(mrl(:,3)>0)'
  for k=1:2
    q{k}(mrl(i,3))=min([q{k}(mrl(i,3)); q_save{k}(mrl(i,1):mrl(i,2))]);
    ads{k}(mrl(i,3))=min([ads{k}(mrl(i,3)); ads_save{k}(mrl(i,1):mrl(i,2))]);
    tmp=[C.dat(mrl(i,3),:); CL21.dat(mrl(i,1):mrl(i,2),:)];
    [v,vi]=max(abs(tmp));
    v2=v.*sign(tmp(vi+(0:size(tmp,1):size(tmp,1)*(size(tmp,2)-1))));
    C.dat(mrl(i,3),:)=v2;
  end
end

regs=regs_save;
fld={'st','en','peak_st','peak_en','peak_wide_st','peak_wide_en','peak','broad_st','broad_en'};
for k=1:2
  for i=1:length(regs{k})
    in_chr=find(C.chrn==CL21_save.chrn(regs{k}(i).peak));
    for j=1:length(fld)
      val=getfield(regs{k}(i),fld{j});
      if ~isempty(val)
        [tmp,nval]=min(abs(pos(in_chr)-CL21_save.pos(val)));
        %      disp([ k i j val tmp nval]);
        regs{k}(i)=setfield(regs{k}(i),fld{j},nval+min(in_chr)-1);
      end
    end
  end
end
