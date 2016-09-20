function history_D(D,show_time)
% history_D displays the recorded history of D
%     history_D(D,show_time)
% Returns no value. If show_time is truthy, the times of history
% records are displayed.
if ~exist('show_time','var')
  show_time=0;
end

if isfield(D,'history')
    for i=1:length(D.history)
        if show_time
            disp([ '(' num2str(i) ')' D.history{i}{1} ': ' D.history{i}{2}]);
        else
            disp([ '(' num2str(i) ')' D.history{i}{2}]);
        end    
        if length(D.history{i})>2
            disp(D.history{i}(3:end))
        end 
    end
end
