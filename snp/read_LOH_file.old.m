function LOH=read_LOH_file(fname)

f=read_dlm_file(fname);

LOH_mat=zeros(length(f)-2,length(f{3})-6);
right_end=length(f{3})-1;
for i=3:size(LOH_mat,1)
  if length(f{i})<2 
    break;
  end
  j=i-2;
  st=cell2mat(f{i}(6:right_end));
  st(st=='L')='1';
  st(st=='R')='2';
  st(st=='N')='0';
  LOH_mat(j,:)=(num2str(st'))';
  LOH_dat(j).marker=str2num(f{i}{1});
  LOH_dat(j).chr=str2num(f{i}{2}); %FIXME should be a string 
  LOH_dat(j).pos=str2num(f{i}{3});
  LOH_dat(j).cM=str2num(f{i}{4});
  LOH_dat(j).score=str2num(f{i}{5});
  if mod(i,1000)==0
    disp(i);
  end
end

LOH_mat=LOH_mat(1:j,:)-48;
LOH_dat=LOH_dat(1:j);

LOH_sdesc=strvcat((f{2}(6:right_end))');

LOH.mat=LOH_mat;
LOH.dat=LOH.dat;
LOH.sdesc=LOH_sdesc;

% save LOH.mat LOH_mat LOH_dat LOH_sdesc
