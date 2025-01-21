function integral = hitomiss(f, m, a, b)
% Calcula l'integral d'una funcio pel metode de montecarlo de hit o miss
% f -referencia a la funcio a la que fer l'integral
% m- nombre de llançaments aleatoris
% a- limit inferior
% b- limit superior
% error: de l'ordre de (1/sqrt(m)
hmax = max(f(a),f(b));
areaRectangle = (b-a)*hmax;
xr = rand(1,m)*(b-a) + a;
yr = rand(1,m)*hmax;
% plot(xr,yr,'o')
favorables = 0;
for ii = 1:m
    if(yr(ii) < f(xr(ii)))
       favorables = favorables + 1;
    end
end

integral = areaRectangle*favorables/m;
end