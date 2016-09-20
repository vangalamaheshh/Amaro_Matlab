function pth=rmpath_regexp(regexp)

p=path;
if ispc
    p1=dlmsep(p,';');
else
    p1=dlmsep(p,':');
end    
p2=grep(regexp,p1);
pth=rmpath(p2{:});
