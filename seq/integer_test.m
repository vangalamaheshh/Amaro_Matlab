boom = '/xchip/tcga_scratch/ng/GBM-0188/wgs/boom/normal.boom';
chr = 2;
st = 1e6;
en = 2e6;
germcall(boom,chr,st,en);




fname = [boom '/chr' num2str(chr) '.base'];


tic,d = get_block(fname,'byte',st,en,'double');toc   % 7.4 sec
tic,b = get_block(fname,'byte',st,en,'uint8');toc    % 5.8 sec
tic,b = get_block(fname,'byte',st,en,'int8');toc     % 5.5 sec

tic,b = get_block(fname,'long',round(st/8),round(en/8),'int64');toc     % 1.3 sec


boom = '/xchip/tcga_scratch/ng/GBM-0188/wgs/boom/normal.boom';
chr = 2;
st = 1e6;
en = 1.5e6;
germcall(boom,chr,st,en);


tic
  offset = B(firstbase,4)-1;
  Kt = sparse(B(firstbase:lastbase,4)-offset,1:(lastbase-firstbase+1),1);
  Lsub = L(:,firstbase:lastbase);
  CL = (Kt*Lsub')';
toc  %Elapsed time is 16.555774 seconds.

  BB = int32(B(:,4));


  offset = B(firstbase,4)-1;
  Kt = sparse(B(firstbase:lastbase,4)-offset,1:(lastbase-firstbase+1),1);
  Lsub = L(:,firstbase:lastbase);
  CL = (Kt*Lsub')';

