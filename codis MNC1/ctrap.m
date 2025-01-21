%Code 6: Composite Trapezoidal Quadrature
% Input: a-b (low-up lim.); m (# intervals); fun (func. name)
% Output: I_{1,m}(f)
function T = ctrap(a,b,m,fun)
    h = (b-a)/m;
    x = a + [0:m]'*h;
    f = feval(fun,x); 
    N = m + 1;
    T = h*(.5*f(1)+sum(f(2:N-1))+.5*f(N));
end