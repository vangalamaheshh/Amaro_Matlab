function gp_scatter_plot(gctfile,rc,dims,clsfile)

if ischar(rc) && ~isempty(str2num(rc))
  rc=str2num(rc);
end

if isnumeric(rc)  
    rc=enum_param(rc,{0,'row';1,'col'});
end

if ischar(dims)
    dims=str2num(dims);
end

R=read_mit_gct_file(gctfile);

if exist('clsfile','var') && ~isempty(clsfile)
    R=read_mit_cls_file(R,clsfile,1);
    c=hsv2rgb([ repmat((0:(1/9):1)',3,1) [ ones(20,1); 0.5*ones(10,1)] ...
        [ones(10,1); 0.5*ones(10,1); ones(10,1);]]);
    c=c(mod((1:7:(size(c,1)*7))-1,size(c,1))+1,:);
    R=add_supmark(R,c);
    R.supmark(1).marker=cellstr(repmat('o',size(R.supmark(1).marker,1),1));
end

if length(dims)<2 || length(dims)>3
    error('Only 2 or 3 dimensions supported');
end

X=R.dat;
if is_col(rc)
    X=X';
end
X=X(:,dims);

figure(1); clf;
if isfield(R,'supdat')
  scatter_plot_D(X,R,1,1);
  [typeacc,typedesc,R1,range,non_empty]=decollapse_supdat(R,1);
  legend(deblank(cellstr(R1.supacc(range,:))));
else
  scatter_plot_D(X,R);
end
grid on

% mcc -m gp_scatter_plot
