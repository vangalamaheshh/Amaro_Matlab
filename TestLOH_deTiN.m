%genreate counts for deTiN LOH
clear
color=jet(8);
power=3;
%for power=1:.5:4
num_points=20;
depth=repmat(100,round(num_points),1);
TiN_V=.25;
AFt=.5;
TiN=[0:0.01:1];

for i=1:num_points

coin=randi([0:1], [1]);
clusters=randi([1,3], [1]);

if clusters==1
    TiN_V=.8;
elseif clusters==2
    TiN_V=.25;
elseif clusters==3
    TiN_V=.1;
end


if coin==1

alt_t(i)=binornd(depth(i),.95);

alt_n(i)=binornd(depth(i),(TiN_V*.5)+.5);

ref_n(i)=depth(i)-alt_n(i);

ref_t(i)=depth(i)-alt_t(i);

else
    
    alt_n(i)=binornd(depth(i),.5-(TiN_V*.5));
    ref_n(i)=depth(i)-alt_n(i);
    
alt_t(i)=binornd(depth(i),.05);

ref_t(i)=depth(i)-alt_t(i);

end

       if coin==1
          pTiN(i,:)=betapdf(.5+(AFt*TiN),alt_n(i)+1,ref_n(i)+1);
          pTiN(i,:)= pTiN(i,:)/sum( pTiN(i,:));
          
       else
           pTiN(i,:)=betapdf(.5-(AFt*TiN),alt_n(i)+1,ref_n(i)+1);
          pTiN(i,:)= pTiN(i,:)/sum( pTiN(i,:)); 
        end
end
pTiN_sum=sum(log(pTiN),1);
pTiN_sum=pTiN_sum+(1-max(pTiN_sum));
plot(exp(pTiN_sum),'Color',color(round(power*2),:))




hold on

pause(.1)
