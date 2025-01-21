function [T,Y] = ab4(f,t0,h,y0,m)
%---------------------------------------------------------------------------
%AB4 Adams-Bashforth solution for y' = f(t,y) with 
% with y(t0) = y0(1) ; y(t0+h) = y0(2).
%   y(t0+2h) = y0(3) ; y(t0+3h) = y0(4).
% Sample call
%   [T,Y] = rk4('f',a,h,ya,m)
% Inputs
%   f    name of the function
%   t0    initial time
%   h    time step
%   y0   initial values
%   m    number of steps
% Return
%   T    solution: vector ( 1 x m )          of abscissas
%   Y    solution: matrix ( length(ya) x m ) of ordinates
%---------------------------------------------------------------------------
  
T = zeros(1,m+1); T(1:4) = t0 + [0 h 2*h 3*h] ;
Y = zeros(length(y0(:,1)),m+1);
Y(:,1) = y0(:,1); F0 = feval(f,T(1),Y(:,1));
Y(:,2) = y0(:,2); F1 = feval(f,T(2),Y(:,2));
Y(:,3) = y0(:,3); F2 = feval(f,T(3),Y(:,3));
Y(:,4) = y0(:,4); F3 = feval(f,T(4),Y(:,4));

for j=4:m
  T(j+1) = t0 + h*j; 
  Y(:,j+1) = Y(:,j) + h*((-3/8)*F0+(37/24)*F1-(59/24)*F2+(55/24)*F3);
  Y(:,j-3) = Y(:,j-2); Y(:,j-2) = Y(:,j-1); Y(:,j-1) = Y(:,j); Y(:,j) = Y(:,j+1);  
  F0 = F1 ; F1 = F2 ; F2 = F3 ; F3 = feval(f,T(j+1),Y(:,j+1));      
end
