%1-ЫЙ ШАГ: ОТБОР СОСЕДЕЙ ДЛЯ КАЖДОЙ ЧАСТИЦЫ
function [x,y] = select(Nb,M,X,Y,Z)

C = 0.75; % радиус отбора соседей для 1 массы

for i = 1:Nb
    k = 1;
    if X(i)==inf && Y(i)==inf && Z(i)==inf
        R(i) = 0;
        J(i,k) = 0;
    else
        R(i) = C+M(i)*(1.025)*(10^(-5));% - радиус отбора соседей
        for j = 1:Nb
            if X(j)==inf && Y(j)==inf && Z(j)==inf
                J(i,k) = 0;  
            else
                if  X(j)-X(i)==0 && Y(j)-Y(i)==0 && Z(j)-Z(i)==0%проверка на то, что j!=i
                    J(i,k) = 0;  
                else
                    if abs(X(j)-X(i))<R(i) && abs(Y(j)-Y(i))<R(i) && abs(Z(j)-Z(i))<R(i)
                        J(i,k) = j; % - вектора номеров соседей для каждой частицы
                        k = k+1;      
                    end
                end
            end
        end
    end
end
x = R;
y = J;
end