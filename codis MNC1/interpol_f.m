function [P,pi,e_n,z] = interpol_f(m,n,fun)
% Input
    % m = nombre de punts de la malla (z)
    % n = nombre de punts equiespaiats (x)
    % fun = funcio a interpolar
% Output
    % P = matriu (m+1) x (n+1) dels polinomis cardinals de Lagrange 
    % evaluats en la malla z. P(i,j) = lambda_i(z_j)
    % pi = vector de m+1 elements resultant de sumar per files els
    % elements de P*f_i, es a dir la interpolació de f en z_j
    % e_n = error maxim comes en cada una de les interpolacions pi
    
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