function write_pvs_regs(fname,C,regs,pvs)

X.marker=C.marker;
X.pos=C.pos;
X.chr=C.chr;
X.chrn=C.chrn;
ampdel={'Amp','Del'};
ad={'A','D'};

X.dat=zeros(length(X.marker),2+4*sum(cellfun('length',regs)));

s=1;
for k=1:length(regs)
  X.sdesc{s}=[ ampdel{k} ': -log10(qv)'];
  X.dat(:,s)=-log10(pvs{k});
  s=s+1;
end


for k=1:length(regs)
  for i=1:length(regs{k})
    % region
    X.sdesc{s}=[ ad{k} 'P' num2str(i) ': region'];
    X.dat(regs{k}(i).st:regs{k}(i).en,s)=1;
    s=s+1;
    
    % wide peak
    X.sdesc{s}=[ ad{k} 'P' num2str(i) ': wide peak'];
    X.dat(regs{k}(i).peak_wide_st:regs{k}(i).peak_wide_en,s)=1;
    s=s+1;

    % peak
    X.sdesc{s}=[ ad{k} 'P' num2str(i) ': peak (qv=' num2str(regs{k}(i).qv) ...
                 ', resid qv=' num2str(regs{k}(i).resid_qv) ')'];
    X.dat(regs{k}(i).peak_st:regs{k}(i).peak_en,s)=1;
    s=s+1;

    % peak center
    X.sdesc{s}=[ ad{k} 'P' num2str(i) ': peak center'];
    X.dat(regs{k}(i).peak,s)=1;
    s=s+1;
    
  end
end

write_as_dchip(fname,X,1);
