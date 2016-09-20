function gp_params=gp_presend_files(gp_params)

global my_gp;
if ~isempty(my_gp)
  f=fieldnames(gp_params);
  for i=1:length(f)
    unix(['mv ' getfield(gp_params,f{i}) ' ' getfield(gp_params,f{i}) '.present']);
    res=runAnalysis(my_gp,'ConvertLineEndings', ...
                    struct('input_filename',[getfield(gp_params,f{i}) '.present'],...
                           'output_file',getfield(gp_params,f{i})));
    gp_params=setfield(gp_params,f{i},res.fileURLs{1});
  end
else
  error('no gp');
end
