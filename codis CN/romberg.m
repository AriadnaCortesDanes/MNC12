% Usage: r = romberg(f,a,b,n)
% Romberg integration
%
% Input:
% f - Matlab inline function 
% a,b - integration interval
% n - number of rows in Romberg tableau
%
% Output:
% r - Romberg tableau containing the computed values of the integral
%
% Examples:
% r=romberg(@sin,0,1,5);
% r=romberg(@myfunc,0,1,5);          here 'myfunc' is any user-defined function in M-file
% r=romberg(inline('sin(x)'),0,1,5);
% r=romberg(inline('sin(x)-cos(x)'),0,1,7);

function r = romberg(f,a,b,n)
h = (b - a) ./ (2.^(0:n-1));
r(1,1) = (b - a) * (f(a) + f(b)) / 2;
for j = 2:n
    subtotal = 0;
    for i = 1:2^(j-2)
        subtotal = subtotal + f(a + (2 * i - 1) * h(j));
    end
    r(j,1) = r(j-1,1) / 2 + h(j) * subtotal;
    for k = 2:j
        r(j,k) = (4^(k-1) * r(j,k-1) - r(j-1,k-1)) / (4^(k-1) - 1);
    end
end
end