%2-ОЙ ШАГ: ВЫЧИСЛЕНИЕ min dt ДЛЯ КАЖДОЙ ЧАСТИЦЫ:
function [r, t]  = search_dt(T,Nb,M,J,X,Y,Z)

kmax = size(J,2);% - max число соседей в массиве
alf = 0.5/T;

for i = 1:Nb
    for j = 1:kmax
        R(i,j) = 0;
    end
end
clear kmax;

for i = 1:Nb
    %в случае если у частицы нет соседей
    if J(i,1) == 0
        minR(i) = 0;% - расстояние между частицей и её ближайшим соседом
        dt(i) = 0;% присваиваем регулярный шаг dT(n,i)
    %в случае если у частицы есть соседи
    else
        k = size(find(J(i,:)),2);% - max число соседей для i-ой частицы
        minR(i) = 0;% - расстояние между частицей и её ближайшим соседом
        %вычисление расстояний между i-ой частицей и её соседями
        for j = 1:k 
            R(i,j) = ((X(J(i,j))-X(i))^2+(Y(J(i,j))-Y(i))^2+(Z(J(i,j))-Z(i))^2)^(1/2);   
            if R(i,j)<minR(i) || minR(i)==0
                minR(i)= R(i,j);
            end
        end
        dt(i) = alf*minR(i)^(3/2);%alf*minR(i)^(3/2);%нашли min иррегурный шаг dt для каждой частицы
        if minR(i)<0.05 && M(i)<1e1
            dt(i) = dt(i)/3;
        end
    end
end
clear reg; 
r = R;
t = dt;
end