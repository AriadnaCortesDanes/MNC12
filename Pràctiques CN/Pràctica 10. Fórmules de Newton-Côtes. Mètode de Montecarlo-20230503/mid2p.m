% Usage: r = mid2p(f,a,b,n)
% Composite 2-point open Newton-Cotes rule
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
% r=mid2p(@sin,0,1,10);
% r=mid2p(@myfunc,0,1,10);          here 'myfunc' is any user-defined function in M-file
% r=mid2p(inline('sin(x)'),0,1,10);
% r=mid2p(inline('sin(x)-cos(x)'),0,1,10);

function r = mid2p(f,a,b,n)
h = (b - a) / (n*3);
hh = h * 2;
x = a + h;
r = 0;
for i=1:n
    r = r + f(x);
    x = x + h;
    r = r + f(x);
    x = x + hh;
end
r = r * h*1.5;
end