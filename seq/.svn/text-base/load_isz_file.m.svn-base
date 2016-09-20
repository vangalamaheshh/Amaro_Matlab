function X = load_isz_file(isz,lanelist,params)
% load_isz_file(isz,lanelist,params)
%
% load and process results of Firehose's InsertSizeByLane and MakeLanelist
%
% Mike Lawrence 2009-10-15

demand_file(isz);

X = [];
tmp = load_matrix(isz);
X.dat = tmp(:,2:end);
X.mx = tmp(end,1);
X.tot = sum(X.dat(:));

if exist('lanelist','var')
  demand_file(lanelist);
  if ~exist('params','var'), params=[]; end
  params = impose_default_value(params,'lane_blacklist','none');

  X.lanes = load_struct(lanelist);

  if strcmpi(params.lane_blacklist,'none')
    BL = [];
  else
    BL = load_lines(params.lane_blacklist);
  end
  X.lanes.is_blacklisted = is_blacklisted(X.lanes.PU,BL);
  
  X = process_isz_data(X,params);
end
