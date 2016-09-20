function Y=match_platforms(Ds,segmented_data)

nDs=length(Ds);

chr_factor=1000;
for i=1:nDs
  chr_factor=max(chr_factor,10.^(ceil(log10(max(Ds{i}.pos)))+2));
end

disp('finding markers');
mrk=[];
n=[];
tot=0;
for i=1:nDs
  % FIXME: add check that each platform is unique and no blanks
  mrk=strvcat(mrk,strvcat(Ds{i}.marker));
  n(i,:)=[ tot+1 tot+size(Ds{i}.dat,1)];
  tot=n(i,2);
  pos{i}=Ds{i}.chrn*chr_factor+Ds{i}.pos;
end
[umrk,umi,umj]=unique(mrk,'rows');
X.marker=cellstr(umrk);

% FIXME: check that pos of same names are the same

posall=cat(1,pos{:});
X.chrn=floor(posall(umi)/chr_factor);
X.pos=mod(posall(umi),chr_factor);

X.dat=sparse(length(umi),1);
X.gsupdat=nan(nDs,length(umi));
for i=1:nDs
  X.gsupdat(i,umj(n(i,1):n(i,2)))=1:(n(i,2)-n(i,1)+1);
end

disp('finding order');
X1=X;
X1.gorigidx=1:size(X1.dat,1);
tmp=zeros(size(X.gsupdat'));
for i=1:nDs
  tmp(:,1:i)=X1.gsupdat(1:i,:)';
  tmp(isnan(tmp))=1e100;
  X1=order_by_pos(X1,tmp);
  
  tmp1.dat=X1.gsupdat(i,:)';
  tmp1=impute_missing_values(tmp1,'before');
  X1.gsupdat(i,:)=tmp1.dat';
end

X2=reorder_D_rows(X,X1.gorigidx);

disp('verify order');
% verify order
for i=1:nDs
  ord{i}=find(~isnan(X2.gsupdat(i,:)));
  tmp=(1:size(Ds{i}.dat,1))-X2.gsupdat(i,ord{i});
  if any(tmp)
    error('Do not match');
  end
end

for i=1:nDs
  Y{i}.marker=X2.marker;
  Y{i}.chrn=X2.chrn;
  Y{i}.pos=X2.pos;
  Y{i}.sdesc=Ds{i}.sdesc;
  if isfield(Ds{i},'sis')
    Y{i}.sis=Ds{i}.sis;
  end
  if isfield(Ds{i},'snp_scores')
    Y{i}.snp_scores=Ds{i}.snp_scores;
  end
  Y{i}.dat=nan(length(X2.chrn),size(Ds{i}.dat,2));
  if exist('segmented_data','var') && segmented_data
    rl=runlength(Ds{i}.dat,Ds{i}.chrn);
    for j=1:length(rl)
      rl{j}(:,1:2)=ord{i}(rl{j}(:,1:2));
    end
    tmp=derunlength(rl);
    Y{i}.dat(1:size(tmp,1),:)=tmp;
  else
    Y{i}.dat(ord{i},:)=Ds{i}.dat;
  end
end


% nS=0;
% crc_table=init_crc;
% for i=1:nDs
% %  nn=crc(uint32(0),crc_table,strvcat(Ds{i}.marker));
% %  nn=str2int_matrix(strvcat(Ds{i}.marker));
% %  nn1=mod(double(nn),1e3);
%   pos{i}=Ds{i}.chrn*chr_factor+Ds{i}.pos; %+nn1/1e3;
%   nS=nS+size(Ds{i}.dat,2);
%   [up,upi,upj]=unique(pos{i});
%   hc=histc(upj,1:length(up));
%   hhc=histc(hc,1:max(hc));
%   exi=find(hc==max(hc));
  
%   if length(up)~=length(pos{i})
%     disp(['Positions are not unique for platform #' num2str(i)]);
%     for j=1:length(hhc)
%       disp([ num2str(j) ': ' num2str(hhc(j))]);
%     end
%     tmp=up(exi(1));
%     disp(['Example: ' num2str(floor(tmp/chr_factor)) ':' sprintf('%ld',mod(tmp,chr_factor))]);
%   end
%   upos{i}=up;
%   uposj{i}=upj;
% end


