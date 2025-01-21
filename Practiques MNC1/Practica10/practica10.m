%% P10 : Cortés y García
clear all
format long

%% Ejercicio 1 : Potencial en el eje
% Función de potencial analítica en el eje
% En y = 0, la integral dada se convierte en:
%       integral de 0 a 2*pi de -d(theta)/sqrt(1+z^2)
%       -(1/sqrt(1+z^2))*(integral de 0 a 2*pi de d(thetha))
%       -2*pi/sqrt(1+z^2)
V_eje = @(w) -2*pi./(sqrt(1+w.^2));

% Calcularemos el potencial en un rango de z = [-5, 5]
z_vec = linspace(-5, 5, 200);

% Para calcular el potencial en (0, 0, z) usamos la función potencial
V_vec = []; V_eje_eval = [];
for z = z_vec
    V_vec = [V_vec potencial(0, 0, z, 1, 1)];
    V_eje_eval = [V_eje_eval V_eje(z)];
end

% Plot de el valor de potencial en z_vec
figure(1)
plot(z_vec, V_vec, 'o', 'Color', 'g')
hold on
plot(z_vec, V_eje_eval, 'Color', 'b', 'LineWidth', 1)
hold off
grid on
legend("Potencial numerico", "Potencial analítico")
title("Potencial en el eje $(0, 0, z)$", 'Interpreter', 'Latex')
xlabel("$z$", 'Interpreter', 'LaTeX')
ylabel("$V$", 'Interpreter', 'LaTeX')

% Observamos que el potencial analitico y el numérico son equivalentes en
% la gráfica, por lo que podemos asumir que la función implementada para el
% cálculo numérico del potencial funciona correctamente.
% Vemos además, como era de esperar, que el comportamiento del potencial en
% el eje z es simétrico respecto a y.

%% Ejercicio 2: Potencial en el plano x= 0

% Repetimos el procedimiento anterior para diferentes valores de y
y_vec = [0.25, 0.75, 1.5];

figure(2)
for y = y_vec
    V_vec = [];
    for z = z_vec
        V_vec = [V_vec potencial(0, y, z, 1, 1)];
    end
    plot(z_vec, V_vec, 'LineWidth', 1);
    hold on
end
hold off
grid on
title("Potencial en $(0, y, z)$", 'Interpreter', 'Latex')
xlabel("$z$", 'Interpreter', 'LaTeX')
ylabel("$V$", 'Interpreter', 'LaTeX')
legend('y = 0.25', 'y = 0.75', 'y = 1.5')

% Si estudiamos el gráfico, sacamos las siguientes conclusiones:
% - Cuando nos encontramos "dentro" del cilindro cuya base viene 
% delimitada por la anilla, el potencial resulta mayor en valor
% absoluto que en puntos externos al mismo.
% - Sin embargo, a medida que el valor de z crece, los potenciales de cada
% caso son más parecidos. Esto se debe a que la anilla pasa a comportarse
% como una masa puntual (para |z| >> 1 o |y| >> 1).

%% Ejercicio 3: Curvas equipotenciales sobre el plano
% Aplicaremos el metodo de diferencias finitas para el calculo del campo
% gravitatorio en diversos puntos iniciales
% El método consiste en que para cada punto, calularemos (Fy, Fz) y lo 
% usaremos para movernos alrededor de la curva equipotencial (usando 
% eps como variación y (Fy, Fz)/norma(Fy, Fz) como dirección) mediante:
% (y1, z1) = (y0, z0) + eps*(Fy, Fz)/norma(Fy, Fz)

% Definimos constantes útiles para el cálculo
h = 0.01; % diferencial que usaremos para las diferencias finitas
eps = 0.001; % diferencial que usaremos para movernos alrededor de las curvas equipotenciales

% Valores iniciales de z que usaremos para dibujar las curvas
% equipotenciales
% Usaremos el valor de y = 0.25 fijado y variaremos z_inicial 
Z_ini = [0:0.1:1]; 

figure(3)
for z_ini = Z_ini
    % Inicializamos los vectores para cada situación inicial
    Y_Vec = [0.25];
    Z_Vec = [z_ini];
    for it = [1:9000]
        y0 = Y_Vec(end);
        z0 = Z_Vec(end);

        Fy = -(potencial(0, y0+h, z0, 1, 1) - potencial(0, y0-h, z0, 1, 1))/(2*h);
        Fz = -(potencial(0, y0, z0+h, 1, 1) - potencial(0, y0, z0-h, 1, 1))/(2*h);
        
        % Aplicamos la iteración dada
        mod = sqrt(Fy^2 + Fz^2);
        y1 = y0 + eps*(-Fz)/mod; 
        z1 = z0 + eps*Fy/mod;

        Y_Vec = [Y_Vec y1];
        Z_Vec = [Z_Vec z1];
    end
    % Representamos la curva de equipotencial
    plot(Y_Vec, Z_Vec, 'LineWidth', 1)
    hold on
end
hold off
grid on
title("Curvas equipotenciales", 'Interpreter', 'Latex')
xlabel("$y$", 'Interpreter', 'LaTeX')
ylabel("$z$", 'Interpreter', 'LaTeX')
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])

% Observamos en la representación gráfica que la curva de nivel de la
% función potencial es simétrica respecto a las dos variables (y, z)

%% Funciones implementadas

% POTENCIAL (potencial.m)
% Input: x,y,z (coordenadas) 
%        lamnda, a (parametros anilla)
% Output: V (potencial en el punto)
% function V = potencial(x, y, z, lambda, a)
%     % Definimos la función del diferencial de potencial
%     dV = @(t) (lambda*a)./sqrt((x-a*cos(t)).^2 + (y-a*sin(t)).^2 + z.^2);
%     
%     V = -1*qclencurt(0, 2*pi, 20, dV);
% end

% QCLENCURT (qclenqurt.m)
% Input: a-b (low-up lim.) 
%        n (# nodes-1) 
%        fun (func. name) 
% Output: I_n(f)
% function Icc = qclencurt(a, b, n, fun)
%     l = [0:n]'; 
%     k = [2:n]'; 
%     x = cos(l*pi/n);
%     
%     w = cos(l*l'*pi/n) \ [2;0;(1+(-1).^k)./(1-k.^2)];
%     
%     z = a+.5*(b-a)*(x+1); 
%     f = feval(fun,z);
%     
%     Icc =.5*(b-a)*w'*f;
% end