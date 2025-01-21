function [P,lamvec] = interpol(m,n)
% Input 
    % m = nombre de punts de la malla (z)
    % n = nombre de punts equiespaiats (x)
% Output
    % P = matriu (m+1) x (n+1) dels polinomis cardinals de Lagrange 
    % evaluats en la malla z. P(i,j) = lambda_i(z_j)
    % lamvec = vector de m+1 elements de la suma dels polinomis l_i per a 
    % cada punt z_j
    
    x = -1 : 2/n : 1; %el vector de punts equiespaiats
    z = -1 : 2/m : 1; %suposem els punts a evaluar son

    P = ones(m+1,n+1); 

    %omplim la matriu amb els coeficients dels productoris
    for ii = 1 : 1 : m+1 %ii representara la fila de la matriu
        for jj = 1 : 1 : n+1 %jj representara la columna de la matriu
             numerador = 1;
             denominador = 1;
             for kk = 1 : 1 : n+1 %kk representara l'index del productori
                 if (jj ~= kk) 
                     numerador = numerador*(z(ii) - x(kk));
                 end
             end
             for kk = 1 : 1: n+1
                 if(jj ~= kk) 
                     denominador = denominador*(x(jj)-x(kk));
                 end
             end
             pi = numerador/denominador;
             P(ii,jj)= pi;
        end  
    end
    lamvec = [];
    
    for ii = 1 : 1 : m+1
        %calcular la funció de Lebesgue per a z(ii)
        lj = P(ii,:);
        lam = sum(abs(lj));
        lamvec = [lamvec lam];
    end
end