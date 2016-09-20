function graphics_helloworld(name,flag)
% graphics_helloworld(name,flag)
%
% testcase for graphics-producing pipeline module
%
% flag=0   --> spots drawn with default color
% flag=1   --> spots drawn multicolored
%
% Mike Lawrence 2009-08-03

if ~exist('name','var'), name = 'world'; end
if ~exist('flag','var'), flag = 0; end
a = rand(200,3);
figure(1);

if flag
  scatter(a(:,1),a(:,2),[],a,'o','filled');
else
  scatter(a(:,1),a(:,2),[],'o','filled');
end

text(0.5,0.5,['Hello ' name '!'],'horizontalalign','center');
print('-dpng','-r600',[name '.png']);
print('-dpdf','-r600',[name '.pdf']);
close(1);
fprintf('Finished successfully.\n');
