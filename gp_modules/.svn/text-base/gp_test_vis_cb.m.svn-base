function gp_test_vis_cb(action)

if ~exist('action','var')
  action='';
end

try
  switch action
   case 'run'
    prev=make_current('A1');
    imagesc(rand(10,10));
    restore_prev('A1',prev);
   otherwise 
    msgbox([mfilename ' ERR: no such action']);
  end
catch
  msgbox([mfilename ' ERR: (' action ') ' lasterr]);
end

function prev=make_current(s)
[co,cf]=gcbo;
prev=gca;
a1_h=findobj(cf,'Tag','A1');
axes(a1_h);

function restore_prev(s,prev)
set(gca,'Tag',s);
axes(prev)
