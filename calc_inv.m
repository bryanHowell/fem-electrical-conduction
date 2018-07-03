function [invA,detA] = calc_inv(A)
% calculate inverse of a square matrix analytically
% -formulas are easily derived and available to public (=> general
% knowledge)
% - analytic formulas are only provided for 3x3 and below

[m,n,o]=size(A);
if(m~=n)
    error('matrix must be square!');
end
if(n>3)
    error('matrix must be 3x3 or smaller!.');
    
end

% pre-allocate
invA = zeros(n,n,o);

if(n==1)
    
    % determinant of matrix
    detA=squeeze(A);

    % inverse of matrix
    invA(n,n,:)=1./A(n,n,:);
    
elseif(n==2)

   % elements of matrix across 3rd dimension 
   a=squeeze(A(1,1,:)); b=squeeze(A(1,2,:));
   c=squeeze(A(2,1,:)); d=squeeze(A(2,2,:));

   % determinant of matrix
   detA=a.*d - b.*c;

   % inverse of matrix
   invA(1,1,:)=d./detA;
   invA(2,2,:)=a./detA;
   invA(1,2,:)=-b./detA;
   invA(2,1,:)=-c./detA;
   
else % DEFAULT
    
   % elements of matrix across 3rd dimension 
   a=squeeze(A(1,1,:)); b=squeeze(A(1,2,:)); c=squeeze(A(1,3,:)); 
   d=squeeze(A(2,1,:)); e=squeeze(A(2,2,:)); f=squeeze(A(2,3,:));  
   g=squeeze(A(3,1,:)); h=squeeze(A(3,2,:)); i=squeeze(A(3,3,:));   
    
   % determinant of matrix.
   detA=a.*(e.*i - f.*h)-b.*(d.*i-f.*g)+c.*(d.*h-e.*g);
   
   % cofactors of matrix inverse
   C11=e.*i - f.*h;
   C12=-(d.*i-f.*g);
   C13=d.*h-e.*g;
   C21=-(b.*i-c.*h);
   C22=a.*i-c.*g;
   C23=-(a.*h-b.*g);
   C31=b.*f-c.*e;
   C32=-(a.*f-d.*c);
   C33=a.*e-b.*d;
   
   % inverse of matrix
   invA(1,1,:) = C11./detA;
   invA(2,1,:) = C12./detA;
   invA(3,1,:) = C13./detA;
   invA(1,2,:) = C21./detA;
   invA(2,2,:) = C22./detA;
   invA(3,2,:) = C23./detA;
   invA(1,3,:) = C31./detA;
   invA(2,3,:) = C32./detA;
   invA(3,3,:) = C33./detA;
   
end

detA=detA(:)'; % consolidate into row vector -> column vector