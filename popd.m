function popd
global SAVED_DIRS
if isempty(SAVED_DIRS)
  disp('No saved dirs.');
else
  cd(SAVED_DIRS{end});
  SAVED_DIRS=SAVED_DIRS(1:(end-1));
end

