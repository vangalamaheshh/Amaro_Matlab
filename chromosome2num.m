function n=chromosome2num(ch)
% n=chromosome2num(ch)
%   replace a chromosome string or array of strings by numbers
%   X = 23, Y = 24, M or MT = 25, XY = 26
%   invalid strings generate NaNs 
%
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.


if isempty(ch)
    n=NaN;
    return;
end
empty_lines=[];
if iscell(ch)
  empty_lines=find(cellfun('isempty',ch));
  if ~isempty(empty_lines)
    ch(empty_lines)=cellstr(repmat('-1',length(empty_lines),1));
  end
  ch=strvcat(ch);
end
ch=lower(ch);
n=zeros(size(ch,1),1);

ch=strvcat(regexprep(cellstr(ch),'chr',''));
%if size(ch,2)>3 && strcmp(ch(1,1:3),'chr')
%  ch=ch(:,4:end);
%end

find_x=strmatch('x',ch,'exact');
if ~isempty(find_x)
  n(find_x)=23;
  ch(find_x,:)=repmat(['0' repmat(' ',1,size(ch,2)-1)],length(find_x),1);
end

find_y=strmatch('y',ch,'exact');
if ~isempty(find_y)
  n(find_y)=24;
  ch(find_y,:)=repmat(['0' repmat(' ',1,size(ch,2)-1)],length(find_y),1);
end

find_m=strmatch('mt',ch,'exact');
if ~isempty(find_m)
  n(find_m)=25;
  ch(find_m,:)=repmat(['0' repmat(' ',1,size(ch,2)-1)],length(find_m),1);
end

find_m=strmatch('m',ch,'exact');
if ~isempty(find_m)
  n(find_m)=25;
  ch(find_m,:)=repmat(['0' repmat(' ',1,size(ch,2)-1)],length(find_m),1);
end

find_m=strmatch('xy',ch,'exact');
if ~isempty(find_m)
  n(find_m)=26;
  ch(find_m,:)=repmat(['0' repmat(' ',1,size(ch,2)-1)],length(find_m),1);
end

% ch=regexp(cellstr(ch),'[^0-9]');
q=double(ch)-48;
q(q~=(32-48) & (q<0|q>9))=NaN;

dig=sum(q~=(32-48),2);
q(q==(32-48))=0;

qq=zeros(size(n,1),1);
for i=1:size(ch,2)
  qq=qq+q(:,i)*10.^(size(ch,2)-i);
end
qq=qq./10.^(size(ch,2)-dig);


%[q,isok]=str2num(ch);
%if ~isok
%   q=NaN*ones(size(ch,1),1);
%   warning('not a chrmosome name');
%end
%n=n+q;

n=n+qq;

if ~isempty(empty_lines)
  n(empty_lines)=NaN;
end
