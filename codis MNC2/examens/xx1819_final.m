%% Final 2018-2019

close all; clear;
format long;

%% Problem 1

%% (b)

mu = 50;

F = @(X) [X(2) ; mu*(1-X(1)^2)*X(2)-X(1)];
Ft = @(t, X) F(X);

X0 = [.01 ; .01];

% Euler Explicit = RK1
hVec = .001:.001:1;
xnormVec = hVec*0;
for k = 1:length(hVec)
    [~, X] = rk1(Ft, 0, hVec(k), X0, 300);
    xnormVec(k) = abs(X(1, end));
end

figure(1)
plot(hVec, xnormVec, '--ok')
% The plot stops ploting when after integrating for T time, the absolute
% value of the final position is infinite (thus NaN).
% The maximum value for F is h=0.011

%% (c)

% Since we want Tm = 100, we'll need m=0000
% [T, X] = rk4(Ft, 0, .0025, X0, 40000);
[T, X] = rk4(Ft, 0, .01, X0, 10000);

figure(2) % Time evolution
plot(T, X(1, :))

figure(3) % Phase space
plot(X(1, :), X(2, :))

%% Problem 2

close all; clear all;
format long;

global lambda alpha beta gamma
global n D D2

%% (a)

fa = @(y) sinh(y)/sinh(1);
yVec = linspace(0, 1, 100);

figure(4)
plot(yVec, fa(yVec));

%% (b)

alpha = 1; beta = -1; gamma = 1; n = 100;
[D, xNodes] = chebdiff(n);

D2 = D^2; yNodes = (1-xNodes)/2;

F0 = @(y) sin(pi*y); f0 = F0(yNodes(2:n));

lambda = 18; [Fk, ~, it] = newtonn(f0, 1e-3, 200, @Lfun);
Fk_18 = Fk(:, it);

lambda = 16; [Fk, ~, it] = newtonn(f0, 1e-3, 200, @Lfun);
Fk_16 = Fk(:, it);

lambda = 14; [Fk, ~, it] = newtonn(f0, 1e-3, 200, @Lfun);
Fk_14 = Fk(:, it);

figure(5)
hold on
plot(yNodes, [0 ; flip(Fk_18) ; gamma ], 'o')
plot(yNodes, [0 ; flip(Fk_16) ; gamma ], 'o')
plot(yNodes, [0 ; flip(Fk_14) ; gamma ], 'o')
hold off

%% (c)

lambdaVec = 18:-1:0;
maxfVec = 0*lambdaVec;
for k = 1:length(lambdaVec)
    lambda = lambdaVec(k); 
    [Fk, ~, it] = newtonn(f0, 1e-3, 200, @Lfun);
    
    maxfVec(k) = max(abs(Fk(:, it)));
end

figure(6)
plot(lambdaVec, maxfVec, 'o');

%% Auxiliar codes

% Problem 1

function [T,Y] = rk1(f, t0, h, y0, m)
    T = zeros(1, m+1);
    Y = zeros(length(y0), m+1);
    T(1) = t0; Y(:,1) = y0;
    
    for j = 1:m
        Y(:,j+1) = Y(:, j) + h*feval(f, T(j), Y(:, j));
        T(j+1) = t0 + h*j;
    end
end

function F = Lfun(f)
    global lambda alpha beta gamma
    global n D D2
    
    c11 = 1; c12 = 0; c13 = 0; 
    c21 = 1; c22 = 0; c23 = gamma;

    L = (4/lambda)*D2 + alpha*eye(n+1);

    M1 = 2*[D(1,2:n) ; D(n+1,2:n)]; 
    M2 = [c21 + c22*(-2)*D(1,1), c22*(-2)*D(1,n+1); c12*(-2)*D(n+1,1), c11+c12*(-2)*D(n+1,n+1)];
    M3 = [L(2:n,1) L(2:n,n+1)];
    
    N = beta*f.*(-2*D(2:n, 2:n)*f); 
    Maux = (M3*inv(M2))*[c22 0; 0 c12]*M1;

    F = (L(2:n, 2:n) + Maux)*f + M3*inv(M2)*[c23 ; c13] + N;
end