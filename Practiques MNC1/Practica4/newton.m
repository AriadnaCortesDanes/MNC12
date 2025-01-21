function [xk,fk,it] = newton(x1,tol,itmax, param)
    it = 0;
    xk = [x1];
    fk = [F(x1, param)];
    ek = 1;
    
    while ek > tol && it < itmax
        x = xk(end) - F(xk(end), param)/derivada(xk(end), param);
        ek = abs(x - xk(length(xk)));
        xk = [xk x];
        fk = [fk F(x, param)];
        it = it + 1;
    end
end