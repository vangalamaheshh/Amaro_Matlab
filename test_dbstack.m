function test_dbstack(n)

[st,i]=dbstack;
disp(st)

if n>1
  test_dbstack(n-1);
end

