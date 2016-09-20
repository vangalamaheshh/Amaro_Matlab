function landscape_plot(varargin)
% Firehose Landscape Plot Module
% get input files
files=options_parsing_for_landscape_plot(varargin);
fields=fieldnames(files);

if isempty(files.sample_to_pair)
error('Sample to Pair mapping is required for plotting')

end


    

end

