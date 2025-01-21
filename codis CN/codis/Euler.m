function [ w,t ] = Euler(f,a,b,h,alpha)
    % referencia: practica 12_estabilitatEuler
    % Mètode d'Euler
    t=[a:h:b];
    N=length(t);
    w(1)=alpha;
    for i=1:N-1
        w(i+1)=w(i)+h*f(t(i),w(i));
    end
end