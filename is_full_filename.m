function res=is_full_filename(fname)

res=0;
if ispc
  if ~isempty(regexp(fname,['(^' regexp_filesep regexp_filesep '|^[a-zA-Z]:' regexp_filesep ')']))
    res=1;
  end
else
  if fname(1)==filesep
    res=1;
  end
end
