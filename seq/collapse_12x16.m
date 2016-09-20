function B=collapse_12x16(A,method)

npages = size(A,3);

if ~(size(A,1)==12 && size(A,2)==16), error('input is wrong size'); end

% first do strand collapse

B=zeros(6,16,npages);
colmap = [16 12 8 4 15 11 7 3 14 10 6 2 13 9 5 1];
for page=1:npages
  for row=1:6
    for col=1:16
      B(row,col,page) = A(row,col,page) + A(13-row,colmap(col),page);
    end
  end
end

A=B;

if strcmpi(method,'washu')
  B=zeros(7,3,2);
  map={[1:4,9:12],[13:16],[5,7,8],[6],[1:12],[13,15,16],[14]};
  for wrow=1:7
    cols=map{wrow};
    for wcol=1:3
      if wrow<=4, row=wcol; else row=wcol+3; end
      for wpage=1:npages
        B(wrow,wcol,wpage) = sum(sum(A(row,cols,wpage)));
      end
    end
  end

elseif strcmpi(method,'broad')
  error('Broad collapse not yet implemented');

elseif strcmpi(method,'strand')
  % no further action needed

else
  error('unknown collapse method');
end

end
