function write_qvs(fname,CL21,pvs,neglog)

Z.marker=CL21.marker;
Z.chr=CL21.chr;
Z.chrn=CL21.chrn;
Z.pos=CL21.pos;
Z.dat=cat(2,pvs{:});
if exist('neglog','var') && neglog
  Z.dat=-log10(Z.dat);
end
Z.sdesc={'Amp','Del'};
write_as_dchip(fname,Z,1,[],1,3);
