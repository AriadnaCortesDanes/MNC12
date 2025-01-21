%Code 7: Composite Simpson’s Quadrature
% Input: a-b (low-up lim.); m (# intervals); fun (func. name)
% Output: I_{2,m}(f)
function S = csimp(a,b,m,fun)
    h = (b-a)/m;
    x = a + [0:2*m]*h/2;
    f = feval(fun,x);
    N = 2*m+1;
    S = (h/6)*(f(1)+4*sum(f(2:2:N-1))+2*sum(f(3:2:N-2))+f(N));
end