function x = load_onecol(infile)
% Mike Lawrence 2009-06-23

t = load_textfile(infile);
x = sscanf(t,'%d');
