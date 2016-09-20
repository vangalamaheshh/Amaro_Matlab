

 D1 =[0    0.5000    1.0300    0.6100; 0.5000         0    0.6700    0.6100;
    1.0300    0.6700         0    0.6200;
    0.6100    0.6100    0.6200         0]
D2 =

         0    0.3536    0.8688    0.5197
    0.3536         0    0.6041    0.4360
    0.8688    0.6041         0    0.4405
    0.5197    0.4360    0.4405         0
Dinf =

         0    0.2500    0.8500    0.5100
    0.2500         0    0.6000    0.3500
    0.8500    0.6000         0    0.3400
    0.5100    0.3500    0.3400         0



function [Val, Ind1, Ind2]=max2d(Mat)
%-------------------------------------------------------------------------- 
% max2d function      2d maximum function. Return the maximum value 
%                   and index in a 2d matrix. 
% Input  : - Matrix. 
% Output : - Maximum value. 
%          - [I,J] index of maximum. 
% Tested : Matlab 5.3 
%     By : Eran O. Ofek                  October 2000 
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html 
%-------------------------------------------------------------------------- 
 
[V1,I1] = max(Mat); 
[V2,I2] = max(V1); 
 
Val = V2; 
Ind1 = I1(I2);   % I 
Ind2 = I2;       % J 
end

function [a, b, Dc] = merge_next( D, c, method )
 
Dc = zeros(length(unique(c)));
 
for i = 1 : size(Dc,1)
    for j = 1 : size(Dc,2)
        Dc(i,j) = calc_distance(i,j,D, c, method);
    end
end
 
Dc( logical(eye(size(Dc,1))) ) = Inf;
[~, a, b] = max2d(-Dc);
 
 
 
 
    function d = calc_distance(i,j, D, c, method)
        switch method
            case 'single'
                d = min(min(D(c == i, c == j)));
            case 'complete'
                d = max(max(D(c == i, c == j)));
            case 'average'
                d = mean(mean(D(c == i, c == j)));
        end
    end
 
 
end
 
 


command line analysis:
Ds = {D1, D2, Dinf}
dists = {'single', 'complete', 'average'}
normnames = {'L1', 'L2', 'Linf'}
 
for i = 1 : length(Ds)
for j = 1 : length(dists)
[a,b] = merge_next(Ds{i}, c, dists{j});
fprintf('Merge: %0.0f with %0.0f\tNorm: %s\tDist: %s\n', a, b, normnames{i}, dists{j});
end
end
 
% OUTPUT
% Merge: 3 with 1   Norm: L1    Dist: single
% Merge: 1 with 1   Norm: L1    Dist: complete
% Merge: 1 with 1   Norm: L1    Dist: average
% Merge: 3 with 1   Norm: L2    Dist: single
% Merge: 1 with 1   Norm: L2    Dist: complete
% Merge: 1 with 1   Norm: L2    Dist: average
% Merge: 3 with 2   Norm: Linf  Dist: single
% Merge: 1 with 1   Norm: Linf  Dist: complete
% Merge: 1 with 1   Norm: Linf  Dist: aver