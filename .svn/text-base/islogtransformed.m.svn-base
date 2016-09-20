function bool = islogtransformed(D,field)
% BOOL = ISLOGTRANSFORMED(D)



chunkdims = getmemchunkdims(D,field,1);

datsize = getsize(D,field);

ii = 0;
jj = 0;

abs_mean = [];

while ii < datsize(2)
  
        thisloopdim2idx = (ii+1):min(ii+chunkdims(2),datsize(2));
        abs_mean(thisloopdim2idx) = mean(abs(D.(field)(:,thisloopdim2idx)),1);
        ii = ii + chunkdims(2);
end
 
    mean_mean = mean(abs_mean);
    
 bool = mean_mean <= .5;