function [c csd] = propagate_error_in_product(a,asd,b,bsd)

c = a.*b;

csd = sqrt(((asd./a).^2)+(bsd./b).^2);
