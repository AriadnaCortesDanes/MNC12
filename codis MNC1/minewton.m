function [xk,fk,it] = minewton(x1,tol,itmax,fun)
% Inputs: 
    % x1 = punt inicial proper a la solucio
    % tol = tolerancia
    % itmax = nombre maxim d'iteracions permeses
    % fun es la funcio de la que volem trobar l'arrel
% Outputs:
    % xk = vector de valor de x que anem trobant fins arribar a l'arrel
    % fk = residu de la funcio fun en cada un dels xk
    
    xk = [x1];
    fk = [fun(x1)];
    ek = 1;
    it = 0;
    
    while ek > tol && it < itmax
        x = xk(end) - fun(xk(end))/derivative(fun,xk(end));
        ek = abs(x - xk(length(xk))); 
        %error comes entre x i l'anterior x trobada
        xk = [xk x];
        fk = [fk fun(x)];
        it = it + 1;
    end
end