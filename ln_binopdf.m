function y=ln_binopdf(x,n,p)

y=chooseln(x,n)+x.*log(p)+(n-x).*log(1-p);  

% [ log(binopdf(0:100,100,0.01)); ln_binopdf(0:100,100,0.01) ]
