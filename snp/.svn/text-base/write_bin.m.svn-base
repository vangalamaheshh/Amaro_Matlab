function write_bin(dirname,C,raise_to_power);
% write binary format
% chr1.array_name.cn.bin
% chr1.snpid.bin
% chr1.snplocation.bin

if(( exist('raise_to_power','var') && raise_to_power) || (~exist('raise_to_power','var') && (abs(mean(C.dat(:)))<0.5)))
  disp('Taking 2.^(x+1)');
  C.dat=2.^(C.dat+1);
end

if exist(dirname,'dir')
  warning('dirname already exists');
else
  mkdir(dirname);
end

dirname=add_slash_if_needed(dirname);
max_chr=max(C.chrn);
for i=1:max_chr
  in_chr=find(C.chrn==i);
  fname=[ dirname 'chr' num2chromosome(i)];
  
  disp(['writing ' fname '.snpid.bin']);
  f=fopen([fname '.snpid.bin'],'w');
  tmp=strvcat(C.marker{in_chr});
  tmp=tmp(:,1:min(size(tmp,2),15));
  tmp=[ tmp repmat(' ',size(tmp,1),15-size(tmp,2))];
  x=zeros(size(tmp,1),30);
  x(:,2:2:(min(size(tmp,2),15)*2))=tmp(:,1:min(size(tmp,2),15));
  fwrite(f,x');
  fclose(f);

  disp(['writing ' fname '.snplocation.bin']);
  f=fopen([ fname '.snplocation.bin'],'w');
  tmp=C.pos(in_chr);
  max_in_chr(i)=max(tmp);
  tmp=uint64(tmp);
  fwrite(f,tmp','uint64',0,'ieee-be');
  fclose(f);

  for j=1:size(C.dat,2)
    disp(['writing ' fname  '.' C.sdesc{j} '.cn.bin']);
    f=fopen([ fname '.' C.sdesc{j} '.cn.bin'],'w');
    tmp=C.dat(in_chr,j);
    tmp=single(tmp);
    fwrite(f,tmp,'single',0,'ieee-be');
    fclose(f);
  end
end

in_chr=1:50:size(C.dat,1);
fname=[ dirname 'chrAll'];
disp(['writing ' fname '.snpid.bin']);
f=fopen([fname '.snpid.bin'],'w');
tmp=strvcat(C.marker{in_chr});
tmp=tmp(:,1:min(size(tmp,2),15));
tmp=[ tmp repmat(' ',size(tmp,1),15-size(tmp,2))];
x=zeros(size(tmp,1),30);
x(:,2:2:(min(size(tmp,2),15)*2))=tmp(:,1:min(size(tmp,2),15));
fwrite(f,x');
fclose(f);

disp(['writing ' fname '.snplocation.bin']);
f=fopen([ fname '.snplocation.bin'],'w');
tmp=double(C.pos(in_chr));
cs=cumsum([0 double(max_in_chr)]);
tmp=(cs(C.chrn(in_chr))')+tmp;
tmp=uint64(tmp);

fwrite(f,tmp','uint64',0,'ieee-be');
fclose(f);

for j=1:size(C.dat,2)
  disp(['writing ' fname  '.' C.sdesc{j} '.cn.bin']);
  f=fopen([ fname '.' C.sdesc{j} '.cn.bin'],'w');
  tmp=C.dat(in_chr,j);
  tmp=single(tmp);
  fwrite(f,tmp,'single',0,'ieee-be');
  fclose(f);
end
