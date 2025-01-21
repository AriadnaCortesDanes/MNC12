%% P5: Cortés y García

clear all
format short
%% Apartado (a) Matriz de polinomios cardinales de Lagrange
n = 10; % numero de puntos equiespaiados 
m = 801; % numero de puntos en los que evaluamos los polinomios
x = -1 : 2/n : 1; %el vector de punts equiespaiats
z = -1 : 2/m : 1; %suposem els punts a evaluar son
 
[A,~] = interpol(m,n);

%% Apartat (b) Casos n = 3, 6, 9
[B1,~] = interpol(m,3);
%plots per columnes (cada un dels polinomis nodals)
figure(1)
for i = 1 : 3
    plot(z,B1(:,i), 'linewidth', 1.5)
    hold on
end
hold off
[B2,~] = interpol(m,6);

figure(2)
for i = 1 : 6
    plot(z,B2(:,i), 'linewidth', 1)
    hold on
end
hold off

[B3,~] = interpol(m,9);
figure(3)
for i = 1 : 9
    plot(z,B3(:,i))
    hold on
end

%% Apartat (c)
figure(4)
[C1,lamvec1] = interpol(m,8);
subplot(2,2,1)
plot(z,lamvec1, 'linewidth', 1, 'color', [0,0.7,0.9])
[C2,lamvec2] = interpol(m,16);
subplot(2,2,2)
plot(z,lamvec2,'linewidth', 1,'color',[0.75,0.25,0.75])
[C3,lamvec3] = interpol(m,24);
subplot(2,2,3)
plot(z,lamvec3,'linewidth', 1,'color',[0.5,0.75,0.5])
[C4,lamvec4] = interpol(m,32);
subplot(2,2,4)
plot(z,lamvec4,'linewidth', 1,'color',[0.93,0.7,0.125])

%% Apartat (d) 
f = @(x) exp(x);

e_maxss = [];
extra = [];
for n = 4 : 2 : 60
    [D, pi, emax] = interpol_f(m,n,f);
    e_maxss = [e_maxss emax];
    extra = [extra, (1e-16)*(2.^n)/(n.*log10(n))]
end

nn = [4:2:60];
error = log10(e_maxss);
figure(5)
plot(nn, error, 'k', 'linewidth', 1)
hold on
plot(nn, log10(extra), 'g', 'linewidth', 1)
hold off





