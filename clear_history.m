function D = clear_history(D)
%CLEAR_HISTORY - remove history fields from D and suppress future history
    if isstruct(D)
        D = rmfield_if_exists(D,{'history','orig','origidx','gorigidx'});
        D.suppress_history = true;
    end
