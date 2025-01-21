%% P7 : Cortes y García

%% Ejercicio 1 - Ejemplo
toeplitz([1 3 5 7],[1 2 3 4])

%% Ejercicio 2,3 
f = @(x) sin(pi.*x);
df = @(x) pi.*cos(pi.*x);

nn = [5 20 40 80];

fig = 1; % Para numerar las figuras

for n = nn
    h = 2/n; 

    % Matriz FD
    v_FD = zeros(1, n+1);
    v_FD(1) = -1;
    v_FD(2) = 1;
    u_FD = zeros(1, n);
    u_FD(1) = -1;
    FD = (1 / h) * toeplitz(u_FD, v_FD);

    % Matriz CD
    v_CD = zeros(1, n+1);
    v_CD(1) = -0.5;
    v_CD(3) = 0.5;
    u_CD = zeros(1, n-1);
    u_CD(1) = -0.5;
    CD = (1/(2 * h)) * toeplitz(u_CD, v_CD);
    
    if n == 5
        disp(FD)
        disp(CD)
    end
    
    % Generamos los nodos equiespaciados
    jj = [0:n]';
    xj = (2 / n) * jj;

    F = f(xj); % evaluación de la función
    dF = df(xj); % evaluación de la derivada de la función

    dF_FD = FD*F; % Aproximación con el método FD
    dF_CD = CD*F; % Aproximación con el método CD

    figure(fig) % Representación para n
    subplot(2, 1, 1)
    plot(xj, dF, 'linewidth', 1);
    hold on
    plot(xj(2:end), dF_FD, 'linewidth', 1, 'color', 'r');
    plot(xj(2:end-1), dF_CD, 'linewidth', 1, 'color', 'g');
    hold off
    grid on
    title('Aproximación de la $df/dx$', 'Interpreter', 'Latex')
    legend('FD', 'CD')

    % Calculamos los errores absolutos respectivos
    e_FD = abs(dF(2:end) - dF_FD);
    e_CD = abs(dF(2:end-1) - dF_CD);

    subplot(2, 1, 2)
    semilogy(xj(2:end), e_FD, 'linewidth', 1, 'color', 'r')
    hold on
    semilogy(xj(2:end-1), e_CD, 'linewidth', 1, 'color', 'g')
    hold off
    grid on
    title('Error de la aproximación de la derivada')
    legend('FD', 'CD')
    
    t = ['Diferenciación numérica para ',num2str(n), ' nodos'];
    sgtitle(t)
    xlim([xj(1), xj(end)])
    
    fig = fig + 1;
end

%% Ejercicio 4

NN = [100:100:2000];

e_max_FD = [];
e_max_CD = [];
hh = [];

for N = NN
    h = 1 / N;
    
    % Matriz FD
    v_FD = zeros(1, N+1);
    v_FD(1) = -1;
    v_FD(2) = 1;
    u_FD = zeros(1, N);
    u_FD(1) = -1;
    FD = (1 / h) * toeplitz(u_FD, v_FD);

    % Matriz CD
    v_CD = zeros(1, N+1);
    v_CD(1) = -0.5;
    v_CD(3) = 0.5;
    u_CD = zeros(1, N-1);
    u_CD(1) = -0.5;
    CD = (1/(2 * h)) * toeplitz(u_CD, v_CD);
    
    % Generamos los nodos equiespaciados
    jj = [0:N]';
    Xj = (2 / N) * jj;

    F = f(Xj); % evaluación de la función
    dF = df(Xj); % evaluación de la derivada de la función

    dF_FD = FD*F; % Aproximación con el método FD
    dF_CD = CD*F; % Aproximación con el método CD
    
    % Calculamos los errores absolutos respectivos
    e_FD = abs(dF(2:end) - dF_FD);
    e_CD = abs(dF(2:end-1) - dF_CD);
    
    % Actualizamos vectores para la representación gráfica
    hh = [hh h];
    e_max_FD = [e_max_FD max(e_FD)];
    e_max_CD = [e_max_CD max(e_CD)];
end

figure(5)
loglog(hh, e_max_FD, 'LineWidth', 1)
hold on
loglog(hh, e_max_CD, 'LineWidth', 1)
hold off
grid on
legend('FD', 'CD')