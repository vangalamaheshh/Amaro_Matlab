function idx=getcoloredtips(R,B)

colored=find((R(:,7)>=2)&(R(:,7)<=15));
intip=false(size(colored));
for j=1:length(colored)
    beg=R(colored(j),9);
    fin=R(colored(j),9)+R(colored(j),5)-R(colored(j),4);    
    intip(j)=(max(sum(B(beg:beg+4,1)>=64), sum(B(fin-4:fin,1)>=64))>2)&(sum(B(beg:fin)==63)<8);
end
idx=colored(intip);