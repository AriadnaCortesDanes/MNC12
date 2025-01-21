function [T,Y] = rk4(f,a,h,ya,m)
%---------------------------------------------------------------------------
%RK4   Runge-Kutta solution for y' = f(t,y) with y(a) = ya.
% Sample call
%   [T,Y] = rk4('f',a,h,ya,m)
% Inputs
%   f    name of the function
%   a    initial time
%   h    time step
%   ya   initial value
%   m    number of steps
% Return
%   T    solution: vector ( 1 x m )          of abscissas
%   Y    solution: matrix ( length(ya) x m ) of ordinates
%---------------------------------------------------------------------------
    T = zeros(1,m+1);
    Y = zeros(length(ya),m+1);
    T(1) = a; Y(:,1) = ya;
    for j=1:m
        tj = T(j); yj = Y(:,j);
        k1 = h*feval(f,tj,yj);
        k2 = h*feval(f,tj+h/2,yj+k1/2);
        k3 = h*feval(f,tj+h/2,yj+k2/2);
        k4 = h*feval(f,tj+h,yj+k3);
        Y(:,j+1) = yj + (k1 + 2*k2 + 2*k3 + k4)/6;
        T(j+1) = a + h*j;
    end