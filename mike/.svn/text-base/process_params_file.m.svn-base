function P = process_params_file(P,fname)

if ischar(P)
  fname = P;
  P = [];
end

if ~isempty(fname)
  if ~exist(fname,'file')
    fprintf('No such params file %s\n',fname);
  else
    tmp=[]; tmp.line = load_lines(fname);
    tmp = parse_in(tmp,'line','^(\S+)\s*(\S*)$',{'key','value'});
%    tmp = load_struct_noheader(fname,2,{'key','value'});
    for i=1:slength(tmp)
      key = tmp.key{i};
      value = tmp.value{i};
      if strcmpi(value,'true'), value=true; end
      if strcmpi(value,'false'), value=false; end
      if ischar(value)
        zzz = str2double(value);
        if ~isnan(zzz), value = zzz; end
      end
      P = impose_default_value(P,key,value);
    end
  end
end
