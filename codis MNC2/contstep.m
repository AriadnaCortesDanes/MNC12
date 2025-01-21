% Code 21: secant continuation step
% Input:    y0 and y1 (two close column vectors)
%           s: pseudo-arclength parameter
%           tol - Newton’s tolerance: ||y_[k+1] - y_[k] || < tol % itmax - max number of iterations
%           fun - function’s name: f(y_1,y_2,...,y_n,y_n+1)
% Output:   y - next point along curve f = 0
%           y belongs to plane orth. to y1-y0
%           passing through secant predictor y1 + s(y1-y0)
%           iconv (0 if y is convergenced to desired tol.)
function [y,iconv] = contstep(fun,y0,y1,s,tol,itmax)
    tolk = 1.0; 
    it = 0; 
    n = length(y0)-1;

    v = y1-y0; 
    yp = y1+s*v; 
    xk = yp;

    while tolk > tol && it < itmax
        Fk = [feval(fun,xk); v'*(xk-yp)]; 
        DFk = [jac(fun,xk); v'];

        [P,L,U] = pplu(DFk); 
        dxk = plusolve(L,U,P,-Fk);

        xk = xk + dxk; 
        tolk = norm(dxk); 
        it = it + 1;
    end

    y = xk;
    if it <= itmax && tolk < tol
        iconv = 0;
    else
        iconv = 1; 
    end
end