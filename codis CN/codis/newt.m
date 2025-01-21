% Usage: r = newt(f,a,b,n)
% Composite Newton's (aka Simpson's 3/8) rule (4-point closed Newton-Cotes)
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
% r=newt(@sin,0,1,10);
% r=newt(@myfunc,0,1,10);          here 'myfunc' is any user-defined function in M-file
% r=newt(inline('sin(x)'),0,1,10);
% r=newt(inline('sin(x)-cos(x)'),0,1,10);

function r = newt(f,a,b,n)
h = (b - a) / (n * 3);
r = f(a);
x = a + h;
for i = 1 : n-1
    r = r + 3 * f(x);
    x = x + h;
    r = r + 3 * f(x);
    x = x + h;
    r = r + 2 * f(x);
    x = x + h;
end
r = r + 3 * f(x);
x = x + h;
r = r + 3 * f(x);
r = r + f(b);
r = r * h*3/8;
end