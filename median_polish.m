function [x,c2,r2,totalEffect]=median_polish(x,max_iter)


if ~exist('max_iter','var')
    max_iter=10;
end

totalEffect=0.0;
delta=0.0;
converged=0;
newSum=0;
oldSum=Inf;
epsilon=0.01;
c=zeros(1,size(x,2));
r=zeros(size(x,1),1);
c2=zeros(1,size(x,2));
r2=zeros(size(x,1),1);
mr=0;
mc=0;

iter=1;
while iter<=max_iter && ~converged
    mc=median(x,1); % median of cols
    x=x-repmat(mc,size(x,1),1);
    c=c+mc;

    % collect totalEffect from c
    delta = median(c);
    c =c - delta;
    totalEffect = totalEffect + delta;
    
    mr=median(x,2); % median of rows
    x=x-repmat(mr,1,size(x,2));
    r=r+mr;

    % collect totalEffect from r
    delta = median(r);
    r =r - delta;
    totalEffect = totalEffect + delta;
    
    iter=iter+1;

    %Check for convergance.
    newSum = sum(sum(abs(x)));
    if (newSum == 0 || abs(1 - oldSum/newSum) < epsilon)
        converged = 1;
%        disp('converged');
        c2 = c + totalEffect;
        r2 = r + totalEffect;
    end;
    oldSum = newSum;

end
return

B = [ 1495.32,1750.07,1809.07,16.904,19.275,9.323,2.978,5.167,4.834,45.501,23.728,17.131;
    1207.07,1248.737,1072.404,16.129,6.898,20.334,3.814,3.831,5.135,28.41,68.57,40.61;
    1054.237,988.154,918.737,47.163,54.988,42.514,13.18,7.341,10.696,11.706,13.035,25.051;
    993.82,850.487,804.154,30.602,30.684,21.785,7.678,8.284,5.489,160.154,93.237,134.737;
    922.987,860.32,974.32,5.361,3.411,5.135,3.444,4.19,4.575,7.935,9.475,14.16;
    855.737,886.487,979.57,7.341,6.542,3.35,5.114,3.647,3.35,13.246,7.831,7.904;
    627.487,586.487,653.904,16.904,4.895,6.818,5.748,8.793,6.612,30.684,7.078,26.731;
    464.654,291.07,350.07,14.495,7.341,18.118,12.433,5.818,7.157,32.95,30.684,45.917;
    453.987,327.404,375.237,3.174,6.542,5.89,2.978,2.835,5.489,11.706,7.831,31.817;
    383.82,430.07,375.237,6.255,6.141,9.735,4.02,7.341,8.221,13.246,11.706,14.16;
    222.487,189.404,214.82,8.221,8.86,10.772,5.404,4.674,4.28,5.748,34.352,2.895];
[X,C,R,T] = median_polish(log2(B),10);
Cres = [9.4171,9.3026,9.3529,3.4967,3.2722,3.6835,2.5231,2.3173,2.6326,4.2899,3.9859,4.4802];
C-Cres

Y=X+repmat(C,size(X,1),1)+repmat(R,1,size(X,2))-T;
median(abs(Y(:)-log2(B(:))))
max(abs(Y(:)-log2(B(:))))


