function m = gkr
% green-black-red colormap

m = zeros(64,3);
m(33:64,2) = 0:(1/31):1;
m(1:32,1) = 1:-(1/31):0;

