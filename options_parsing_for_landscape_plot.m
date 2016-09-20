function [files]=options_parsing_for_landscape_plot(varargin)

% gathers files for plot construction
% options must be input in this way:
% 'Type of File 'path to file'
% Possible file types include:
% ['sample_to_pair' two column file with pair_ids and sample_id (required!!)
% 'Abs_maf' 'Germ_maf' 'absolute_table' 'sig_genes' 
% 'gistic_table' 'sample_order']
% amaro@broadinstitute


for i=1:length(varargin{1})
   if iseven(i)
       if isequal(varargin{1}{i},'-')
           disp(strcat('setting_',current_input,'_to empty'))
        files.(current_input)=[];
       else
           disp(strcat('Loading_',current_input))
           files.(current_input)=load_struct(varargin{1}{i});
       end  
   else
   current_input=varargin{1}{i};
   end
   
    
    
end






end