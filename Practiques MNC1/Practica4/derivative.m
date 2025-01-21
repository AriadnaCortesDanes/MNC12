function der = derivative(f, x0)
    dx = 1e-6;
    der = (f(x0+dx)-f(x0))/dx;
end