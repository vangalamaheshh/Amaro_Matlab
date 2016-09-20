function X = get_library_list(basedir)
% get_library_list(basedir)
%
% traverses the specified directory and looks for all files matching the form:
% Solexa-*.hybrid_selection_metrics
%
% returns a struct with the following fields:
%    library = Solexa-*
%    sample = TCGA-*
%    version = v*
%    baitset = taken from *.hybrid_selection_metrics column "BAIT_SET"
%
% Mike Lawrence 2009-07-14

X = [];
X.library = {};
X.sample = {};
X.version = {};
X.baitset = {};

X = get_library_list_recurse(X,basedir);

  function X = get_library_list_recurse(X,basedir);
%    fprintf('Checking %s\n',basedir);
    d = dir(basedir);
    for i = 1:length(d)
      if d(i).name(1)=='.', continue; end
      fname = [basedir '/' d(i).name];
      if d(i).isdir
        X = get_library_list_recurse(X,fname);
      else
        if ~isempty(regexp(d(i).name,'^Solexa-(\d*)\.hybrid_selection_metrics$','tokens'))
          x = parse({fname},'^.*(TCGA-.*)/(v.*)/(Solexa-.*)/Solexa-.*\.hybrid_selection_metrics$',...
            {'sample','version','library'});
          txt = load_lines(fname);
          x.baitset = {regexprep(txt{end},'^(\S*).*$','$1')};
          X = concat_structs({X,x});
        end
      end
    end
  end
end


