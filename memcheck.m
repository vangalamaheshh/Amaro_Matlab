% function memcheck
if(1)
  tmps=whos;
  for tmpi=1:length(tmps);
    if tmps(tmpi).bytes>1e7
      disp([tmps(tmpi).name ':' num2str(tmps(tmpi).bytes)]);
    end
  end
else
  whos
end

unix('free -g');
