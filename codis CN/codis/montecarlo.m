function integral = montecarlo(f, m, a, b)
% Calcula l'integral d'una funcio pel metode de montecarlo del valor
% esperat
% f -referencia a la funcio a la que fer l'integral
% m- nombre de llançaments aleatoris
% a- limit inferior
% b- limit superior
% error: de l'ordre de (1/sqrt(m)
r = rand(1,m)*(b-a) + a;
integral = (b-a)*sum(f(r))/m;
end