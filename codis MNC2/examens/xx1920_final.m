%% Final 2019-2020

close all; clear;
format long;

%% Theoretical exercise

%% (b)

rhoVec = [1 ; -300/137 ; 300/137 ; -200/137 ; 75/137 ; -12/137];
roots_rho = roots(rhoVec);

% disp(roots_rho)
% disp(norm(roots_rho(1)))
% disp(norm(roots_rho(2))) % same as disp(norm(roots_rho(3)))
% disp(norm(roots_rho(4))) % same as disp(norm(roots_rho(5)))

%% (c)

rho = @(mu) mu.^5 + (-300*mu.^4 + 300*mu.^3 - 200*mu.^2 + 75*mu - 12)/137;
sigma = @(mu) (mu.^5)*60/137;
z_ratio = @(mu) rho(mu)./sigma(mu);

N = 150; 
thVec = 2*pi*(0:N-1)/N;
muVec = exp(1i*thVec');
zVec = z_ratio(muVec); % compute the marginal values of z

A = [-1 -7 ; 7 -1];
eVal = eig(A); % disp(eVal);
hVec = .001:0.001:2;

figure(1)
plot(real(zVec), imag(zVec))
hold on
plot(real(hVec*eVal(1)), imag(hVec*eVal(1)), '--k')
plot(real(hVec*eVal(2)), imag(hVec*eVal(2)), '--k')
hold off

% After seeing the lines in the picture, we will want to use Newton's
% method to find the actual intersection values, or approximate them. To
% use Newton's method, we would need to implement a function that
% represents de distance between the line and the boundary.
% The integrations will be numerically stable out of the region.

%% Problem 1

%% (a)

% Let X = [r, th]
% X_t = F(X)

[eVec1, eVal1] = eig(jac(@ F1, [0 ; 0])); eVal1 = diag(eVal1);
[eVec2, eVal2] = eig(jac(@ F1, [1 ; 0])); eVal2 = diag(eVal2);
[eVec3, eVal3] = eig(jac(@ F1, [.5 ; -.5])); eVal3 = diag(eVal3);

disp(eVal1); disp(eVal2); disp(eVal3);
% The second fixed point has real eigenvalues
disp(eVec2);

%% (b)

Ft = @(t, X) F(X);

figure(2)
plot([0 ; 1 ; .5], [0 ; 0 ; -.5], '*k')

hold on
for y0 = -1:0.1:1
    X0 = [1.5 ; y0];
    [~, X] = rk4(@ Ft1, 0, 5e-2, X0, 1000);

    plot(X(1, :), X(2, :));
end

for y0 = -.4:-0.03:-.49
    X0 = [.5 ; y0];
    [~, X] = rk4(@ Ft1, 0, 5e-2, X0, 1000);

    plot(X(1, :), X(2, :));
end
hold off

%% (c)

Y0 = -.4:-0.01:-.49;
Tm = Y0*0; m = 200;

figure(3)
plot(.5, -.5, '*k')
hold on
for k = 1:length(Y0)
    X0 = [.5 ; Y0(k)];
    [T, X] = rk4(@Ft1, 0, 5e-2, X0, m);
    
    dist = zeros(1, m+1);
    for ii = 1:m+1
        dist(ii) = norm(X(:, ii)-X0);
    end
    [~, idx] = min(dist(1:end)); Tm(k) = T(idx+9);

    plot(X(1, :), X(2, :));
end
hold off

disp(Tm)

%% (d)

% No estamos usando eigenvectors :(

[T, trajVec] = rk4(@Ft1, 0, 5e-2, [1, -1e-6], 600);

figure(3)
plot(1, 0, '*k')
hold on
plot(trajVec(1, :), trajVec(2, :))
hold off

%% Problem 2

close all; clear;
format long;

%% (a)

global n D D2 P R
n = 100; 
[D, xNodes] = chebdiff(n);
sNodes = (1 - xNodes)/2; D2 = D^2;

%% (b)

Th0 = sNodes(2:n)*0;

P = 10; R = 5;
[Thk1, ~, ~] = newtonn(Th0, 1e-6, 200, @F2);

R = 10;
[Thk2, ~, ~] = newtonn(Th0, 1e-6, 200, @F2);

%% (c)

figure(3)
hold on
plot(sNodes(1:n), [0 ; flip(Thk1(:, end))])
plot(sNodes(1:n), [0 ; flip(Thk2(:, end))])
hold off

%% (d)

grid = linspace(0, 1, 100);
thetaVec = interp1(sNodes(2:n), Thk1(:, end), grid);

figure(4)
plot(grid, thetaVec)

xVec = thetaVec*0; yVec = thetaVec*0;
for kk = 1:99
    xVec(kk) = trapz(cos(thetaVec([1 kk])));
    yVec(kk) = trapz(sin(thetaVec([1 kk])));
end

figure(5)
plot(xVec, yVec, 'o')

%% Auxiliar codes

function F = F1(X)
    r2 = @(X) X(1)^2 + X(2)^2;
    Fun = @(X) [-X(1)*r2(X) - X(2)*(1-r2(X))+X(1)^2-X(2)^2 ; ...
        X(1)*(1-r2(X)) - X(2)*r2(X) + 2*X(1)*X(2)];
    
    F = Fun(X);
end

function F = Ft1(t, X)
    F = F1(X);
end

function F = F2(f)
    global n D2 D P R
    
    c11 = 1; c12 = 0; c13 = 0; 
    c21 = 0; c22 = 1; c23 = 0;

    L = 4*D2; Lh = L(2:n, 2:n);

    M1 = 2*[D(1,2:n) ; D(n+1,2:n)]; 
    M2 = [c21 + c22*(-2)*D(1,1), c22*(-2)*D(1,n+1); c12*(-2)*D(n+1,1), c11+c12*(-2)*D(n+1,n+1)];
    M3 = [L(2:n,1) L(2:n,n+1)] ; 
    
    N = P*sin(f) - R*cos(f); 
    Maux = (M3*inv(M2))*[c22 0; 0 c12]*M1;

    F = (Lh + Maux)*f + M3*inv(M2)*[c23 ; c13] + N;
end