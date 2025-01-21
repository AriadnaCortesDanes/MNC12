%% Pràctical 8

clear; close all;
format long G;

%% (a) The ODE
% If we let $u(y, t)=G(y)F(t)$, when we compute the EDO resulting from the
% previous equation, we have:
%
% $\frac{d^2F}{dt^2}G = g\frac{d}{dy}\{yF\frac{dG}{dy}\}
% = gF\{\frac{dG}{dt} + y\frac{d^2G}{dt^2}\}$
%
% This can be written as:
%
% $\frac{1}{F} \frac{d^2F}{dt^2} = \frac{g}{G} \{ y\frac{d^2G}{dy^2} + \frac{dG}{dy}\}$
%
% This will only have a solution if both sides of the equation are equal to
% a constant $k$. Now, depending on the sign of $k$, we'll have different
% solutions. Now if we take the temporal component of the wave function we
% have $\frac{1}{F(t)} \frac{d^2F(t)}{dt^2} = k$
%
% - For $k < 0$, let $k = -\lambda ^2$ (with $\lambda \in R$). We have that
% $F(t) = Ae^{i\lambda t}+Be^{-i\lambda t}$. 
% 
% - For $k = 0$, we have $F(t) = Ct + D$.
% 
% - For $k > 0$, let $k = \lambda ^2$ (with $\lambda \in R$). We then have
% that $F(t) = Fe^{\lambda t}+Ge^{\lambda t}$.
%
% The former two options ($k \geq 0$) are not bounded $\forall t$, since
% for $t \to \infty$ both solutions tend to $\infty$ as well. Therefore, 
% they will not be suitable solutions to our problem. However, the 
% solution for $k < 0$ is feasible. Then, let us rewrite the resulting 
% equation for $G(y)$ by using $k=-\lambda^2$.
%
% $y\frac{d^2G}{dy^2} + \frac{dG}{dy} + \frac{\lambda^2}{g}G = 0$
%
% which is exactly the equation we were looking for.

%% (b) Eigenvalues & Eigenfunctions
% From the statement of the practical we can define a the boundary 
% condition for $y=l$: $u(l,t) = 0$, which means $G(l) = 0$. However, even
% though we know that $u(y, t)$ must be bounded for $y=0$, we have no
% defined value a concrete boundary condition for $y=0$ (it is the free end
% of the chain).

% Define the parameters of the problem
n = 26;
[D,x] = chebdiff(n); D2  = D*D; % 1st and 2nd differentiation matrix
y = (x+1)/2; % map the chebyshev nodes to the actual domain of y

%%
% Now, from the equation given in (a), we have that 
% $y\frac{d^2G}{dy^2} + \frac{dG}{dy} = -\frac{\lambda^2}{g}G$. Let us
% write $L = y\frac{d^2}{dy^2} + \frac{d}{dy}$ and 
% $\mu = -\frac{\lambda^2}{g}$.
%
% Since there is only one boundary condition for $G$ and it is $G(l)=0$, we
% would have $c_{11}=1$ and $c_{12}=c_{13}=c_{21}=c_{22}=c_{23}=0$.
% Therefore, it doesn't really make sense to compute the operators $M_1$,
% $M_2$ and $M_3$ since they will be of no use. Now, we'll want to solve 
% $LG = \mu G$ by finding the eigenvalues $\mu$ and eigenvectors $G$ of 
% the ODE. 

L = 4*diag(y)*D2 + 2*D; % Define L operator
L = L(2:end, 2:end); % Take out the column corresponding to y=l

% Find eigenvectors and eigenvalues of the operator
[V, E] = eig(L); E = diag(E);
[E, Perm] = sort(E); 
V = V(:, Perm);
lambda = sqrt(-E); lambdafirst = lambda(end-2:end);

disp(lambdafirst); % Display the smallest 3 eigenvectors

%%
% Plot the three eigenfunctions of the smallest eigenvalues
figure(1)
subplot(1, 3, 1)
hold on
plot([0; V(:, end)], y, 'color', [0.5 0 0])
xline(0, '--k');
hold off
grid on
xlim([-0.5 0.1])
title('$G_3(y)$ for $\lambda_1=1.204$ ', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter','latex')

subplot(1, 3, 2)
hold on
plot([0; V(:, end-1)], y, 'color', [0 0.5 0])
xline(0, '--k');
hold off
grid on
xlim([-0.3 0.5])
title('$G_2(y)$ for $\lambda_2=2.76$ ', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter','latex')

subplot(1, 3, 3)
hold on
plot([0; V(:, end-2)], y, 'color', [0 0 0.5])
xline(0, '--k');
hold off
grid on
xlim([-0.3 0.5])
title('$G_3(y)$for $\lambda_3=4.326$ ', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter','latex')

sgtitle('Eigenfunctions of smallest negative $\lambda$s','Interpreter', 'latex')

%% (c) Bessel functions & nodes
% The nodes of the wave are those values of $y$ for which $u(y, t)=0$. We
% also know that $J_0(2\lambda_k)=0$. Then, we may be able to find the
% nodes of each oscillation mode by looking for 
% $2\lambda_i \sqrt(y)=2\lambda_k$. This means we will be able to compute
% the nodes for the oscillation mode $k$ with 
% $y=\left(\frac{\lambda_i}{\lambda_k}\right)^2$ for $i \leq k$. We'll
% plot this positions in their respective plots of $G(y)$.

for kk = 1:3
    for ii = 1:kk
        zero = (lambdafirst(kk)/lambdafirst(ii))^2; % compute zero
        zero_round = round(zero,2); % just for graphical purposes
        
        figure(1)
        subplot(1, 3, (4-ii))
        if abs(zero -1) > eps
            text(0,zero,['\quad $y =$ ' num2str(zero_round)], 'Interpreter', 'latex')
        end
        hold on
        plot(0, zero, 'ok');
        hold off
    end
end

%% Auxiliar codes

% Code 5B: Chebyshev Differentiation matrix
% Input: n 
% Output: differentiation matrix D and Chebyshev nodes
function [D,x] = chebdiff(n)
    x =cos([0:n]'*pi/n); d = [.5 ;ones(n-1,1);.5];
    D = zeros(n+1,n+1); 
    for ii = 0:n
     for jj = 0:n
      ir = ii + 1 ; jc = jj + 1;
      if ii == jj
       kk = [0:ii-1 ii+1:n]'; num = (-1).^kk.*d(kk+1) ;
       D(ir,jc) =((-1)^(ir)/d(ir))*sum(num./(x(ir)-x(kk+1)));
      else
       D(ir,jc) = d(jc)*(-1)^(ii+jj)/((x(ir)-x(jc))*d(ir));
      end
     end
    end
end