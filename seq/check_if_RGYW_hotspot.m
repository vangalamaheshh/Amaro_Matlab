function h = check_if_RGYW_hotspot(v)
% h = check_if_RGYW_hotspot(v)
%
% v is 5bp "vicinity" reference sequence centered on mutated base

h = false(length(v),1);
for i=1:length(v)
  q=upper(v{i});
  if q(3)=='C'
    if q(1)=='A' || q(1)=='T'
      if q(2)=='G' || q(2)=='A'
        if q(4)=='C' || q(4)=='T'
          h(i)=true;
    end,end,end
  elseif q(3)=='G'
    if q(2)=='G' || q(2)=='A'
      if q(4)=='C' || q(4)=='T'
        if q(5)=='A' || q(5)=='T'
          h(i)=true;
    end,end,end
  end
end
