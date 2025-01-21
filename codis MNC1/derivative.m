function [puntprima] = derivative(fun, punt)
% Numerical derivative of fun at punt
    delta = 10e-8;
    f0 = fun(punt);
    f1 = fun(punt + delta);
    puntprima = (f1 - f0)/delta;
end


