function host=get_hostname

[s,r]=unix('hostname');

if s==0
    host=r;
else
    host='Error';
end 
