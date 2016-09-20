function r = geometric_series(st,en,num)

f = en/st;
step = f.^(1/(num-1));

r = as_column(st*step.^(0:num-1));
r(1) = st;
r(end) = en;
