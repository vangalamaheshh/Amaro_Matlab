%Tetrad analysis
%e.g. take phenotype table :
% 0:4 1:3 2:2 3:1 4:0 and output scenario :)
% 0 0 15 66 19
% wt x mutant



% y1 = input('What is the first organism e.g. wt mut \n' ,'s'); % right now set up for only WT x Mut
% y2 = input('What is the second organism \n','s');
y1='MUT'; y2='WT';
display(sprintf('The Cross is %s with %s \n',y1,y2));

%type_of_analysis = input('What are we testing e.g. tetrad , punnet , etc. \n','s'); %implementing tetrad first.
type_of_analysis='tetrad';
if strfind(type_of_analysis,'tetrad')>0
    disp('Doing tetrad analysis')
end
% Get inputs

display('Please input the phenotype ratios from left 0:4 (all wt) to right 4:0 (all mut) \n')

counts_matrix=input('');

PD=counts_matrix(3); % the parental type is 2:2
TT=counts_matrix(2)+counts_matrix(4); % the tetratype is 3:1 or 1:3
NPD=counts_matrix(1)+counts_matrix(5); % the NPDs are all 0:4 or 4:0
N=sum(counts_matrix);

if binocdf(PD,N,NPD/N)>.05 && NPD>0 % check if PD equals NPD ; binomial sorting alleles.
    disp('PD equals NPD genes are unlinked')
    if TT > 0 % there are tetratypes
        D=((.5*TT)/N)*100;
        disp(sprintf('If A is known to be tightly linked then this is the distance of B: %i',round(D)));
        
        X2=((TT-(4/6*N))/(4/6*N)+(PD-(1/6*N))/(1/6*N)+(NPD-(1/6*N))/(1/6*N)); % X2 table for TT,PD,NPD given by A being 100% linked and B being as unlinked as possible
        if chi2cdf(X2,2,'upper') > .05
            disp('This is likely a 1:1:4 ratio meaning that one gene is 100% centromere linked and the other is very far away')
            return
        else
            disp('Looks like this is 1:1:<4 ratio this indicates that at least one gene is not 100% centromere linked')
            return
        end
    else
        disp('There are no tetratypes therefore this both are 100% linked to different centromeres')
        return
    end
    
    
else
    disp('PD does not equal NPD ')
    if TT > 0 % there are tetratypes
        disp('PD does not equal NPD and there are tetratypes therefore these genes are linked')
        D=((.5*(6*NPD+TT))/N)*100;
        disp(sprintf('Distance in cM is %d ',D));
        return
    elseif TT == 0 % redundant for clairity
        if NPD > 0
            disp('PD does not equal NPD and there are no tetratypes therefore these genes are 100% linked distance = 0cM')
            return
        else
            disp('All tetrads are PD therefore this is one gene that is sorting in a Mendelian fashion')
            return
        end
    end
    
end

