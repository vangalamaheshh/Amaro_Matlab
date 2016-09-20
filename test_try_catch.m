close all
try
  w=get(1);
catch
  disp(1);
  try
    e=get(1);
  catch
    disp(2);
  end
end
