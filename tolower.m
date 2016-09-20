function lst=tolower(st)

lst=st;
uppos=find( (st >='A') & (st <='Z'));
lst(uppos)=char(double(st(uppos)))-(double('A')-double('a'));


