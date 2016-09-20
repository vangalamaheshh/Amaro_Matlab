function s=write_line(txt,dlm)

s=[];
if iscell(txt)
  s=txt{1};
  if ~ischar(s)
    s=num2str(s);
  end
  for i=2:length(txt)
    si=txt{i};
    if ~ischar(si)
      si=num2str(si);
    end
    s=[ s dlm si ];
  end
else
  s=deblank(txt(1,:));
  for i=2:size(txt,1)
    s=[ s dlm deblank(txt(i,:))];
  end
end  
