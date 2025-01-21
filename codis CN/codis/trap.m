% Usage: r = trap(f,a,b,n)
% Composite trapezoid rule (2-point closed Newton-Cotes)
%
% Input:
% f - Matlab inline function 
% a,b - integration interval
% n - number of subintervals (panels)
%
% Output:
% r - computed value of the integral
%
% Examples:
% r=trap(@sin,0,1,10);
% r=trap(@myfunc,0,1,10);          here 'myfunc' is any user-defined function in M-file
% r=trap(inline('sin(x)'),0,1,10);
% r=trap(inline('sin(x)-cos(x)'),0,1,10);

function r = trap(f,a,b,n)
h = (b - a) / n;
r = f(a) * 0.5;
x = a + h;
for i = 1 : n-1
    r = r + f(x);
    x = x + h;
end
r = r + f(b) * 0.5;
r = r * h;
end