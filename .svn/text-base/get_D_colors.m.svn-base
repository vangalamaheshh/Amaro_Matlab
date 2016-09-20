function cols=get_D_colors(D,rc,supid)

if is_col(rc)
    cols=zeros(size(D.dat,2),3);
    [u,ui,uj]=unique_keepord(D.supdat(supid,:));
    cm=D.supmark(supid).colormap;
    if size(cm,1)<max(uj)
        error('too few colors');
    else
        cols=cm(uj,:);
    end
else
    cols=zeros(size(D.dat,1),3);
    [u,ui,uj]=unique_keepord(D.gsupdat(supid,:));
    cm=D.gsupmark(supid).colormap;
    if size(cm,1)<max(uj)
        error('too few colors');
    else
        cols=cm(uj,:);
    end
end