%% Codis robats del llibre que no tinc clar que fan maybe

% solve a equation of the form 
%       p(x)f'' + q(x)f' + r(x)f = g(x)
% with bounded domain conditions 
%       c11 f(-1) + c12 f'(-1) = c13
%       c21 f(1) + c22 f'(1) = c23
% using Chevyshev nodes to discretize the system of equations

n = 32 ;
[D,x] = chebdiff(n) ;
D2 = D*D ; 
I = eye(n+1);
Q = diag(-pi*cos(pi*x)) ;
L = D2 + Q*D + I ;
c11 = 1; c12 = -1; c13 = 1+pi; c21 = 1; c22 = 1; c23 = 1-pi;
g = exp(sin(pi*x)).*(1-pi^2*sin(pi*x));

M1 = -[D(1,2:n) ; D(n+1,2:n)] ; 
M3 = [L(2:n,1) L(2:n,n+1)];
M2 = [c21 + c22*D(1,1), c22*D(1,n+1);  c12*D(n+1,1), c11+c12*D(n+1,n+1)];
N = L(2:n,2:n) + (M3*inv(M2))*[c22 0;0 c12]*M1 ;
f = N\(g(2:n)-M3*inv(M2)*[c23 ; c13]);

plot(x(2:end-1),f,'o')
title('Numerical aproximation at Chebyshev nodes with n = 32')

%% Eigenvalue problem
% this is the same idea but for solving the particular case
%       p(x)f'' + q(x)f' + r(x)f = lambda·f
% with bounded domain conditions 
%       c11 f(-1) + c12 f'(-1) = 0
%       c21 f(1) + c22 f'(1) = 0
% using Chevyshev nodes to discretize the system of equations

n = 26 ; 
[D,x] = chebdiff(n) ; 
D2 = D*D ; 
I = ones(n+1);
c11 = 1 ; c12 = 0 ; c13 = 0 ; c21 = 1 ; c22 = 2 ; c23 = 0 ;

L = 4*D2 ; M1 = -[D(1,2:n) ; D(n+1,2:n)] ;
M2 = [c21 + c22*D(1,1), c22*D(1,n+1);
c12*D(n+1,1), c11+c12*D(n+1,n+1)];
M3 = [L(2:n,1) L(2:n,n+1)] ;
N = L(2:n,2:n) + (M3*inv(M2))*[c22 0;0 c12]*M1 ;

[EVEC,EVAL] = eig(-N) ; lamb = diag(EVAL) ;
[foo,ii] = sort(lamb) ; lamb = lamb(ii) ; EVEC = EVEC(:,ii) ;

plot(x(2:end-1),EVEC(:,1)), hold on,
plot(x(2:end-1),EVEC(:,2))
plot(x(2:end-1),EVEC(:,3))
legend('$f_1(z)$', '$f_2(z)$', '$f_3(z)$', 'Interpreter', 'latex')
title('Eigenfunctions associated with $\lambda_1, \lambda_2, \lambda_3$','Interpreter', 'latex')


%% Same shit but non-linear 
%   L f = N(f,g)    % N holds for the non-linear part of the equation
% with initial function f_0(x) = 1 
% using chebyshev differentiation
% mirar codi de la Laura si surt aixo :) -> en serio

global n D D2 I
n = 32 ; 
[D,x] = chebdiff(n) ; 
D2 = D*D ; 
I = ones(n+1);
itmax = 24 ; 
tol = 1e-12 ;
f0 = ones(n-1,1);
xn = x(2:n);
[XK,resd,it] = newtonn(f0,tol,itmax,@fsb);
y =.5*(x+1);
plot(y(2:n),XK(:,end),'-ok'); hold on
title('Solution of the non-lineal BVP amb condicio incial f0')

% funcio que calcula les matrius linials i no lineals i esas cosas
% exemple: f'' = 1.5 f^2
% function F = fsb(f)
%     global n D D2 I
%     L = D2; c11 = 1; c12 = 0; c13 = 4; c21 = 1; c22 = 0; c23 = 1;
%     M1 = -[D(1,2:n) ;D(n+1,2:n)] ; 
%     M3 = [L(2:n,1) L(2:n,n+1)] ;
%     M2 = [c21 + c22*D(1,1), c22*D(1,n+1);
%     c12*D(n+1,1) , c11+c12*D(n+1,n+1)];
% 
%     N = (3/8)*f.^2;
%     Maux = (M3*inv(M2))*[c22 0;0 c12]*M1;
%     
%     F = (L(2:n,2:n) + Maux)*f + M3*inv(M2)*[c23;c13] - N;
% end

%% Same shit per periodic domains
% -f'' + 4 cos(2x) f = E f
clear all;

N = 64 ; x = 2*pi*[0:N-1]'/N; D2 = dftdiffmat2(N);
q = 4.0*cos(2*x) ; L = -D2 + diag(q) ;
[EVEC,EVAL] = eig(L) ; 

lamb = diag(EVAL) ;
[foo,ii] = sort(lamb) ; lamb = lamb(ii) ; EVEC = EVEC(:,ii) ;

plot(x(:,1),EVEC(:,1)), hold on,
plot(x(:,1),EVEC(:,2))
plot(x(:,1),EVEC(:,3))
legend('$f_1(z)$', '$f_2(z)$', '$f_3(z)$', 'Interpreter', 'latex')
title('Eigenfunctions associated with $E_1, \E_2, \E_3$','Interpreter', 'latex')

%% Cosas de unbounded domain
n = 64 ; [D,x] = chebdiff(n) ; D2 = D*D ;
L0 = 5 ; Q = 1-x.^2; V = x.^2;
L = (-diag(Q.^3)*D2+diag(3*x.*Q.^2)*D)/L0^2+diag(L0^2*V./Q);
N = L(2:n,2:n) ; [EVEC,EVAL] = eig(N) ;
