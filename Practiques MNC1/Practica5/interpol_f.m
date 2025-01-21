function [P,pi,e_n] = interpol_f(m,n,fun)
    x = -1 : 2/n : 1;
    z = -1 : 2/m : 1;
    
    [P,~] = interpol(m,n); % calculamos P con la función que tenemos
    
    fx = fun(x); % evaluamos la función en los polinomios modales
    pi =[]; % polinomio interpolador
    
    for i = 1 : m+1
        pi_act = sum(fx.*P(i,:));
        pi = [pi pi_act];
    end
    
    fz = fun(z);
    error = abs(pi - fz);
    e_n = max(error);
end