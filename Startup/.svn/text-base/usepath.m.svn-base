function usepath(varargin)
%USEPATH  create new user source root(s) for matlab path.
% (remove all previous paths down to $HOME/matlab, keep matlab system paths)
% Takes multiple string arguments.

% only modify matlab path if not in MCR environment
if ~isdeployed

    % validate arguments before changing matlab path
    for i = 1:length(varargin)
        root = varargin{i};
        if ~ischar(root)
            error(['USEPATH arguments must be strings; argument ',num2str(i),' is not.']);
        end
        if ~exist(root,'dir')
            error('matlab:usepath:notadir','''%s'' is not a directory!',root);
        end   
    end
    
    % remove path down to $HOME/matlab (startup.m path) 
    fprintf('Removing all user paths\n');
    x = regexp(path,':','split');
    hi = find(strcmp([getenv('HOME') '/matlab'],x),1);
    if length(hi)==1 && hi ~= 1
        rmpath(x{1:hi-1});
    end
    % add each path and subdirectories to tha matlab path
    for i = 1:length(varargin)
        root = varargin{i};
        fprintf('Adding directories under %s to matlab path...\n',root);
        addtree(root);
    end
end
