function riker_make_freeze(S,params)
% riker_make_freeze(S,params)
%

rikerdir = '/xchip/cga1/lawrence/riker';

if ~exist('S','var'), error('requires "S"'); end

if ~exist('params','var'), params = []; end
params = impose_default_value(params,'WGS_only',false);
params = impose_default_value(params,'capture_only',false);
params = impose_default_value(params,'freeze_file','freeze.txt');

if params.WGS_only & params.capture_only, error('Incompatible parameters'); end

if params.WGS_only
  S = reorder_struct(S,S.tech==1);
elseif params.capture_only
  S = reorder_struct(S,S.tech==2);
end

fprintf('Saving freeze as %s\n',params.freeze_file);
save_struct(S,[rikerdir '/' params.freeze_file]);
