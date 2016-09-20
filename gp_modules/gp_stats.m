function gp_stats(infile,outfile,rc)

if ischar(rc) && ~isempty(str2num(rc))
  rc=str2num(rc);
end

if isnumeric(rc)  
    rc=enum_param(rc,{0,'row';1,'col'});
end

D=read_mit_gct_file(infile);
S=calc_D_stats(D,rc);

write_mit_gct_file(outfile,S);

% mcc -m gp_stats

