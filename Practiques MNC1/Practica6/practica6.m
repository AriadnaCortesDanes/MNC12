%% P6: Cortes y García
clear all
format long

%% Apartado (a) Representación de la función
f = @(x) tanh(20*sin(12*x)) + (2/100).*exp(3*x).*sin(300*x);

xx = linspace(0,1,2000);
figure(1)
plot(xx, f(xx),'linewidth',1,'color',[0 0.7 0.9])
title('Representacion de f(x)','Interpreter','Latex','fontSize',16)
xlabel('$x$','Interpreter','Latex')
ylabel('$f(x)$','Interpreter','Latex')

%% Apartado (b) Interpolación con distinto numeros de nodos
figure(2)
e_max = []; %vector que almacena el error maximo cometido con cada valor de n
for n = 100:50:2000
    j = [1:n-1]; % no incluimos j = 0 y j = n el vector para tratarlos aparte
    z = cos((pi * j) / n); % Nodos de Chebyshev
    x = 0.5*(1 + z); % Cambio de variable
    % en j = 0 -> z = 1 -> x = 1
    % en j = n -> z = -1 -> x = 0
    
    interp = []; % vector de interpolaciones
    
    %interpolacion por la forma baricéntrica
    for x_val = xx
        % Tratamos los casos de los nodos j = 0 (x = 0) y j = n (x = 1)
        sum_num = 0.5.*f(1)/(x_val-1) + 0.5.*f(0)/(x_val);
        sum_den = 0.5/(x_val-1) + 0.5/(x_val);
        
        % Tratamos el resto de nodos
        for jj = j
            sum_num = sum_num + (-1)^jj.*f(x(jj))/(x_val - x(jj));
            sum_den = sum_den + (-1)^jj/(x_val - x(jj));
        end
        
        interp = [interp sum_num/sum_den];
    end
    error = abs(interp - f(xx));
    e_n = max(error);
    e_max = [e_max e_n];
    
    % hacemos el plot solo de n = 100, 200, 300...
    if mod(n,100) == 0
        e = log10(error);
        plot(xx, e);
        hold on
    end
end
hold off
title('Error cometido en $[0,1]$ para distintas n','Interpreter','Latex','fontSize',18)
xlabel('$x$','Interpreter','Latex','fontSize',16)
ylabel('$\epsilon(x)$','Interpreter','Latex','fontSize',16)
legend('$n=100$', '$n=200$', '$n=300$', '$n=400$','$n=500$', '$n=600$', '$n=700$', '$n=800$', '$n=900$', '$n=1000$', '$n=1100$', '$n=1200$','$n=1300$', '$n=1400$', '$n=1500$', '$n=1600$', '$n=1700$', '$n=1800$', '$n=1900$', '$n=2000$','Interpreter','Latex')

% Podemos apreciar que el error local cometido por la interpolación
% disminuye a medida que n aumenta. Además, observamos picos en los valores
% de x en los que, al representar f(x), vemos que esta varia su punto 
% central de oscilación drásticamente entre 1 y -1. Esto tiene sentido,
% puesto que una pequeña variación de x en uno de estos puntos supone un
% gran cambio en el valor de f(x), lo cual lleva a un mayor error local.

%% Apartado (c) Numero mínimo de nodos para un error concreto
% Aprovechamos el bucle anterior para calcular el error local máximo para
% cada n, de manera que solo es necesario buscar dicho valor
error_tol = 1e-6; % error máximo 
n_min = find(e_max <= error_tol);
100 + 50*(n_min(1)-1) % n_min(1) es un indice asi que encontramos el valor que le corresponde