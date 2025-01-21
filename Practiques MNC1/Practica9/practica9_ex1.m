%% P9: Cortés y Gracía
clear all
format long

%% Ejercicio 1: Exactitud
f1 = @(x) 8*x.^7;
f2 = @(x) 5*x.^4;
f3 = @(x) 4*x.^3;
a = 0; b = 1; % valors dels extrems de la integral
error1 = []; %vectores que guardaran els errors comesos per cada valor de m
error2 = [];
error3 = [];
for mm = 1: 10 :80 
    H = (b-a)/mm;
    k = [0:1:2*mm];
    xk = a + k*H/2;
    
    %calculem els diferents valors que necesitem per a cada integral
    f01 = f1(xk(1));
    sum11 = sum(f1(xk(3:2:2*mm-1)));
    sum21 = sum(f1(xk(2:2:2*mm)));
    f2m1 = f1(xk(2*mm+1));
    
    I1 = H/6*(f01+2*sum11+4*sum21+f2m1);
    error1 = [error1 abs(1-I1)];

    %per a la segona integral 
    f02 = f2(xk(1));
    sum12 = sum(f2(xk(3:2:2*mm-1)));
    sum22 = sum(f2(xk(2:2:2*mm)));
    f2m2 = f2(xk(2*mm+1));
    
    I2 = H/6*(f02+2*sum12+4*sum22+f2m2);
    error2 = [error2 abs(1-I2)];
    
    %per a la tercera integral 
    f03 = f3(xk(1));
    sum13 = sum(f3(xk(3:2:2*mm-1)));
    sum23 = sum(f3(xk(2:2:2*mm)));
    f2m3 = f3(xk(2*mm+1));
    
    I3 = H/6*(f03+2*sum13+4*sum23+f2m3);
    error3 = [error3 abs(1-I3)];
end

%fem un plot en escala logaritmica dels errors comesos en funcio de m a
%cada integral
figure(1)
loglog([1: 10 :80], error1, 'lineWidth', 1)
hold on
loglog([1: 10 :80], error2, 'lineWidth', 1)
hold on
loglog([1: 10 :80], error3, 'lineWidth', 1)
hold off
grid on
legend('f = 8x^7','f = 5x^4','f = 4x^3')

%vemos que el error decae a medida que augmentamos el valor de m.
%Observamos también que el error es menor cuando el grado de f es menor.
%Remarcar que en el caso de f de grado menor 4 el error tiene directamente
%a 0 (0 numeric de matlan), ya que como vemos analiticamente el error es proporcional 
%a la derivada cuarta de la funcion, que en estos casos es nula. En los casos de
%polinomios de grado >= 4, el error decrece cuando crece m (teoricamente como m^(-4)).
%comprovamos que efectivamente el error decae como m^(-4) con la funcion polyfit
polyfit(log10([1: 10 :80]), log10(error1),1)
polyfit(log10([1: 10 :80]), log10(error2),1)
polyfit(log10([1: 10 :80]), log10(error3),1)

%% Ejercicio 2: Aproximando pi

%% Ejercicio 3: Opcional