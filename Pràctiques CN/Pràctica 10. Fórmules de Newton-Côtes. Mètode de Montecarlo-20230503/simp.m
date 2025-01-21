% Usage: r = simp(f,a,b,n)
% Composite Simpson's rule (3-point closed Newton-Cotes)
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
% r=simp(@sin,0,1,10);
% r=simp(@myfunc,0,1,10);          here 'myfunc' is any user-defined function in M-file
% r=simp(inline('sin(x)'),0,1,10);
% r=simp(inline('sin(x)-cos(x)'),0,1,10);

function r = simp(f,a,b,n)
h = (b - a) / (n * 2);
r = f(a);
x = a + h;
for i = 1 : n-1
    r = r + 4 * f(x);
    x = x + h;
    r = r + 2 * f(x);
    x = x + h;
end
r = r + 4 * f(x);
r = r + f(b);
r = r * h/3;
end