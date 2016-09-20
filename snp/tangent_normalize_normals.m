function [Y]=tangent_normalize_normals(A)  %#ok
% TANGENT_NORMALIZE_NORMALS  Apply tangent normalization to all columns of
% a matrix using leave-one-out for each column of the matrix.
%
%       Y = TANGENT_NORMALIZE_NORMALS(A)
%
%       Each column k in Y is the minimized distance (in the least squares sense)
%       between the kth column of A and the closest vector in the range of
%       the rest of A.  
%
%       Uses qrinsert and qrdelete to allow for computation of Q,R only
%       once. (Question: Does this mean that samples normalized toward the
%       end of the iteration will have accumulated more floating point
%       error?)
%
%           Revisions:
%               8 May 08 -- Function added.  (jdobson@broad.mit.edu)

%---
%$Id$
%$Date: 2008-07-29 11:45:35 -0400 (Tue, 29 Jul 2008) $
%$LastChangedBy: jdobson $
%$Rev$




%% TEST IF THIS WORKS BY ADDING A LINEARLY DEPENDENT ROW TO INPUT DATA
%% FIRST, rotate the normals to a space that is square

[Qreduce,Rreduce] = qr(A,0);
Areduce = inv(Rreduce)*(pinv(Qreduce)*A);

%% NOW normalize using qrinsert/qrdelete

% Leave out first column, then cycle through columns leaving out ith column
% on ith iteration

N = Areduce;


% QR decomposition
disp('QR DECOMPOSITION')
tic
[Q,R]=qr(N,0);
toc
disp('DONE WITH QR')
%Determine r = effective rank 
tol = eps*norm(N,'fro');
r = sum(abs(diag(R)) > tol);
r = r(1); % fix for case where R is vector.
% % Use first r columns of Q.
% if r > 0
%     Q = Q(:,1:r);
%     R = R(1:r,:);
% else
%     Q = [];
% end

Y = zeros(size(Areduce));

% Now normalize the rest using qrdelete and qrinsert
for k = 1:size(Areduce,2)
    
    disp(k)
    
    disp('QR DELETE')
    [Q1,R1] = qrdelete(Q,R,k);


    disp('NORMALIZING')
    Y(:,k) = Areduce(:,k) - Q*(Q'*Areduce(:,k));
    
    if mod(k,100)
        disp(num2str(k))
        toc
    end
end

    
fprintf(1,'\n')
    
Y = Qreduce*Rreduce*Y;
% 
% 
% % Leave out first column, then cycle through columns leaving out ith column
% % on ith iteration
% 
% N = A(:,2:end);
% 
% 
% % QR decomposition
% disp('QR DECOMPOSITION')
% tic
% [Q,R]=qr(N,0);
% toc
% disp('DONE WITH QR')
% %Determine r = effective rank 
% tol = eps*norm(N,'fro');
% r = sum(abs(diag(R)) > tol);
% r = r(1); % fix for case where R is vector.
% % % Use first r columns of Q.
% if r > 0
%     Q = Q(:,1:r);
%     R = R(1:r,:);
% else
%     Q = [];
% end
% 
% keyboard
% disp('Done inverting R.  Normalizing first normal');
% Y = zeros(size(A));
% Y(:,1) = A(:,1) - Q*(Q'*A(:,1));
% 
% 
% % Now normalize the rest using qrdelete and qrinsert
% for k = 2:size(A,2)
%     
%     disp(k)
%     
%     disp('QR DELETE')
%     Q = qrdelete(Q,R,k-1);
%     disp('QR INSERT')
%     Q = qrinsert(Q,R,k-1,A(:,k-1),'col');
% 
%     disp('NORMALIZING')
%     Y(:,k) = A(:,k) - Q*(Q*A(:,k));
%     
%     if mod(k,100)
%         disp(num2str(k))
%         toc
%     end
% end
% 
%     
% fprintf(1,'\n')
%     

 
 
 
 