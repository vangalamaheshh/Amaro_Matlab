function C = distinct_colors(n)
% distinct_colors(n)
%
% Mike Lawrence 2009-07-15

if isnan(n)   % special legacy case
  C = [1 0.6 0.6;0.6 1 0.6;0.6 0.6 1;1 1 0.6;1 0.6 1;0.6 1 1;1 1 1];
  return
end

Q = [...
0 0 0;... % black
1 0 0; 0 0 1; 0 0.6 0;...  % red dk.blue dk.green
1 0.5 0; 0.6 0 0.6; 0 0.7 0.9;... % orange purple lt.blue
0.9 0.8 0; 1 0.7 0.7; 0.4 1 0.4;...   % yellow pink lt.green
0.7 0.7 0.7;...   % grey
];

nq = size(Q,1);

if n<=nq, C = Q(1:n,:);
else
  C = rand(n,3);
  C(1:nq,:) = Q;
end

