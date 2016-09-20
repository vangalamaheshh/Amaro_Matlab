function quick_disp(X)

if size(X.dat,1)>1e5
  imagesc_trim(X.dat(1:10:end,:));
else
  imagesc_trim(X.dat);
end

bluepink;
