 TiN=0:.01:1;
 
for i=1:slength(call)
    pTiN(i,:)=betapdf(tumor_f*TiN),n_alt_count+1,n_ref_count(i)+1);
end

TiN_bin = 1./( 1 + 1*((TiN/tumor_f) - 1) )

TiN_bin=2./((b.tau(seg)./TiN)-b.tau(seg)+2);
pTiN=histw(pTiN,TiN_bin,TiN);