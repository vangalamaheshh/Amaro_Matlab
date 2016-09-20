function D=fix_D(D)


if ischar(D.gdesc)
  D.gdesc=strs2cell(D.gdesc);
end

if ischar(D.gacc)
  D.gacc=strs2cell(D.gacc);
end


