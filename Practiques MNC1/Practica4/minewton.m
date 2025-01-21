function [xk,fk,it] = minewton(x1,tol,itmax,fun)
    it = 0;
    xk = [x1];
    fk = [fun(x1)];
    ek = 1;
    
    while ek > tol && it < itmax
        x = xk(end) - fun(xk(end))/derivative(fun,xk(end));
        ek = abs(x - xk(end));
        xk = [xk x];
        fk = [fk fun(x)];
        it = it + 1;
    end
end