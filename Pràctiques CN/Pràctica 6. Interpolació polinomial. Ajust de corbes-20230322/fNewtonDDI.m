function [ L,sol ] = fNewtonDDI( x,y,xf )
%   x absicises
%   y dades a interpolar
%   punts on evaluar el polinomi interpolador
n=length(x); % Number of terms of X or Y
L=zeros(n);% For divided difference functions.
for k=1:n
   L(k,1)=y(k);
end
for k=2:n 
  for m=1:n-k+1                     
       L(m,k)=(L(m,k-1)-L(m+1,k-1))/(x(m)-x(m+(k-1)));
  end
end
%disp(L);
s=1;t=0;
sol=y(1);
for k=2:n
   s=s*(xf-x(k-1));
   t=s*L(1,k);
   sol=sol+t;
end
end