function guess=heuristic_min2(ctr,Z,testset,trainset,Q,debug)

% calc Field
if nargin==5
  debug=0;
end
Ztest=Z(trainset,testset);
ntest=size(Ztest,2);
Field=sparse(ntest,Q);% matrix of all candidates verses the optional classification
for i=1:Q
  typei=find(ctr(trainset)==i);
  if length(typei)==1
    Field(:,i)=Ztest(typei,:)'; %only one-> no need to sum
  else
    Field(:,i)=sum(Ztest(typei,:))';
  end
end
TotField=sum(Field'); % sum of field on each point
posF=find(TotField>0); % list of points with a positive field


J=Z(testset,testset); %all the neigbours
clear Z
% avJ=mean(nonzeros(J));
guess=zeros(ntest,1);

UnknownField=sum(J); % sum of the J to unknown neighbors
noF=find(TotField==0); % list of point with no field
nt=ntest; % nt - # left to be classified
alone=intersect(find(UnknownField==0),noF); % point with no neighbors = alone
                                            % deal with points with no hope (no neighbors)
if length(alone) > 0
  guess(alone)=-1;
  UnknownField(alone)=-1;
  nt=ntest-length(alone);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trying to classify all unclassified points %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while( nt > 0 ) % go over all unclassified points
  
  if debug 
    disp('begining');
    [ 1 nt ntest ]
  end
  [dum1,S]=max(Field'); % find for each point the largest field 
  
  % classify all points with no unknown neighbors
  unclass=find(guess==0);
  % small number because of roundoff errors in UnknownField
  noUnknowns=unclass(find(UnknownField(unclass)<=0.001)); % noUnknowns-all the points with ratio->infinity (obvious classification)
  guess(noUnknowns)=S(noUnknowns)'; % classify them to the type with largest field
  UnknownField(noUnknowns)=-1;
  nt=nt-length(noUnknowns);
  
  if ( nt > 0 ) %if there are more unclassified left
    
    % calculate delta
    F=Field';
    F(S+(0:ntest-1)*Q)=0; %hide the maximal -> can find the second best
    [dum2,S2]=max(F);
    delta=dum1-dum2;
    
    % calc ratio=delta/unknown field and choose max
    rat=delta./UnknownField;
    [dum,site]=max(rat);
    if debug
      [ 2 nt ] 
      dum
    end
    if dum > 1 % if there is at least 1 with rat > 1
      sites=find(rat>1); 
      guess(sites)=S(sites); % classify all points with rat>1 to the maximal field
      UnknownField(sites)=-1;
      nt=nt-length(sites);
      % update the field and unknown-field of its neighbors
      [dum,nb,j]=find(J(sites,:));
      if length(nb) > 0
        for ii=1:length(nb)
          Field(nb(ii),S(sites(dum(ii))))=Field(nb(ii),S(sites(dum(ii))))+j(ii);
          UnknownField(nb(ii))=UnknownField(nb(ii))-j(ii);		       
        end
      else
        warning('should not reach this point 1')
        keyboard
      end
      
    else % there is no point with rat > 1. Perform the heuristic step
         %		[nt length(find(guess==0)) dum]
      if dum == 0 %the maximal ratio is zero: two cases  1.no neibours 2. the equel rats-confused
                  % confused between several equal options -> choose first
        confused=intersect(find(dum1~=0), ...
                           find(guess==0));
        if length(confused)>0        
          guess(confused)=S(confused);
          %			    [5 confused S(confused) dum1(confused) dum2(confused)]
          %			    Field(confused,:)
          UnknownField(confused)=-1;
          nt=nt-length(confused);
          % no clue -> reject
          %  nt=nt-length(find(guess==0));
          
          [dum3,nb,j]=find(J(confused,:));
          if length(nb) > 0
            for ii=1:length(nb)
              Field(nb(ii),S(confused(dum3(ii))))=Field(nb(ii),S(confused(dum3(ii))))+j(ii);
              UnknownField(nb(ii))=UnknownField(nb(ii))-j(ii);		       
            end
          else
            warning('should not reach this point 2')
            keyboard
          end
        else
          nt=0;
          guess(find(guess==0))=-1;    
        end
      else % dum > 0 => we have a prefered candidate
           % classify maximal
        if guess(site) ~= 0 
          disp('BAD 2')
          keyboard
        end
        guess(site)=S(site);
        UnknownField(site)=-1;
        nt=nt-1;
        [dum3,nb,j]=find(J(site,:));
        if length(nb) > 0
          Field(nb,S(site))=Field(nb,S(site))+j';
          UnknownField(nb)=UnknownField(nb)-j;
        else
          warning('should not reach this point 3');
          keyboard
        end
      end
    end
  end	
end







