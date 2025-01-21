% Code 1: Bisection method for solving f(x) = 0 in [a,b]
% Input: 1. [a,b]: interval (it assumes that f(a)f(b) < 0)
% 2. tol: tolerance so that abs(x_k+1 - x_k) < tol
% 3. itmax: maximum number of iterations allowed
% 4. fun: function’s name
% Output: 1. xk: resulting sequence
% 2. res: resulting residuals
% 3. it: number of required iterations
function [xk,res,it] = bisection(a,b,tol,itmax,fun)
    ak = a; 
    bk = b; 
    xk = [];
    res = [];
    it = 0; 
    tolk = abs(bk-ak)/2 ;
while it < itmax & tolk > tol
    ck = (ak + bk)/2;
    xk = [xk ck];
    if it > 0; tolk = abs(xk(end)-xk(end-1)); end
    fa = feval(fun,ak); fc = feval(fun,ck); res = [res abs(fc)];
    if fc*fa < 0; bk = ck; else ak = ck; end
    it = it + 1 ;
end