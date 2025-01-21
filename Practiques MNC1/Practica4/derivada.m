function der = derivada(x0, param);
    dx = 1e-6;
    der = (F(x0+dx, param)-F(x0, param))/dx;
end