%   Entrada: valor de x, vector con parámetros incluidos en la definición
%   de la función
%   Salida: valor de la función evaluada en x dados los parámetros

function [f] = F(x, param)
    %param es un vector tal que param[1] es v0 y param[2] es alpha0
    g = -9.81;
    v0 = param(1);
    alpha0 = param(2);
    a = param(3);
    b = param(4);
    
    f = x.*tand(alpha0)+0.5*g*(x/(v0*cosd(alpha0))).^2-a.*x.^2.*exp(-b.*x);
end