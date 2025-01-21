% Usage: r = midp(f,a,b,n)
% Composite midpoint rule (1-point open Newton-Cotes)
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
% r=midp(@sin,0,1,10);
% r=midp(@myfunc,0,1,10);          here 'myfunc' is any user-defined function in M-file
% r=midp(inline('sin(x)'),0,1,10);
% r=midp(inline('sin(x)-cos(x)'),0,1,10);

function r = midp(f,a,b,n)
h = (b - a) / n;
x = a + h * 0.5;
r = 0;
for i=1:n
    r = r + f(x);
    x = x + h;
end
r = r * h;
end