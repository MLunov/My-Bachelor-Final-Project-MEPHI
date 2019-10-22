%3-ИЙ ШАГ: ПОИСК ЧАСТИЦ С САМЫМ МАЛЕНЬКИМ min dt(иррегулярным временным шагом)
function [iMIN,t,S,cluster,out_claster,bwn,bwon]  = mindt(dt,M,J,num_b)

s=1;%счетчик частиц с min dt для которых ведется перерасчет
dtmin(s) = min(dt(find(dt~=0)));%dtmin - самый min dt;  
k = find(dt==dtmin(s));
p1(s) = k(1);%p1 - номер частицы с минимальным dt
% [dtmin(s),p1(s)] = min(dt);%dtmin - самый min dt; p1 - номер частицы с минимальным dt 

t_int = dt; % - intermediate dt - промежуточный вектор шагов dt
%отбор частиц с наименьшим min dt, фиксируем их соседей и dt
while find(t_int)~=0 
    %Вычисление окончательного №  частицы для которой будет расчёт      
    %считаем суммарную массу всех соседей частицы 
    k = find(J(p1(s),:));
    k1 = size(k,2);
    Mclust1(s)=M(p1(s));
    p2(s)=0;
    for i = 1:k1
        Mclust1(s)=Mclust1(s)+M(J(p1(s),i));
        %находим частицу-пару у которой dtmin тот же
        if t_int(J(p1(s),i)) == dtmin(s)
            p2(s) = J(p1(s),i);
        end
    end
    if p2(s) == 0
        Mclust2(s)=0;
    else
        %считаем суммарную массу соседей частицы-пары i
        k = find(J(p2(s),:));
        k2 = size(k,2);
        Mclust2(s)=M(p2(s));
        for i=1:k2 
            Mclust2(s)=Mclust2(s)+M(J(p2(s),i));
        end
    end
    %сравниваем массы их кластеров
    %получаем окончательно №  частицы для которой будет расчёт
    if Mclust2(s)>Mclust1(s)
        imin(s) = p2(s);
        k = k2;
        Miclust(s) = Mclust2(s);
    else
        imin(s) = p1(s);
        k = k1;
        Miclust(s) = Mclust1(s);
    end
    clear k1; clear k2;

    % удаление из массива № соседей i
    for i = 1:k 
        Cl(s,i)=J(imin(s),i);%-зафиксировали вектора номеров соседей каждой частицы с min dt
        nb = J(imin(s),i);%-номер частицы соседа
        t_int(nb)=0;%-зануляем их временной шаг из-за их ненадобности
        J(nb,:)=0;%-обнуляем у соседей вектора номеров соседей
        %обнуляем номера соседей в других векторах
        [q,m] = find(J==nb); 
        a=find(q~=imin(s));
        l = find(a);
        l = size(l,1);
        for j=1:l
            u = find(J(q(a(j)),:));
            u = size(u,2);
            J(q(a(j)),m(a(j)))=0;
            for v = 1:(u-m(a(j)))
                J(q(a(j)),m(a(j))-1+v)=J(q(a(j)),m(a(j))+v);
            end
            J(q(a(j)),u)=0;
            if  J(q(a(j)),1)==0% - условие отсутствия соседей
                t_int(q(a(j)))=0;
            end
        end
    end
    clear nb; clear q;  clear m;  clear a;  clear l;  clear u;  clear v; 

    t_int(imin(s))=0;
    J(imin(s),:)=0;
    

    %Снова находим частицу с минимальным иррегулярным шагом dt
    %поиск первого dtmin
    k = find(t_int); 
    if k~=0
        s = s + 1;
        l = size(k,2);
        dtmin(s)= t_int(k(1));
        p1(s) = k(1);

        for i=1:l
            if  t_int(k(i))<dtmin(s)
                dtmin(s) = t_int(k(i));
                p1(s) = k(i);%номер частицы с наименьшим минимальным dt
            end
        end
    end
end
clear t_int;  
iMIN=imin;
t=dtmin;
% cl = Cl;
S = s;
%добавление номера imin-ой частицы к вектору номеров её соседей
CL = Cl;
for z = 1:s
    i = imin(z);%переобозначим номер частицы с наименьшим dt
    k = find(Cl(z,:));
    k = size(k,2);%-число соседей
    CL(z,k+1) = i;%вектора номеров частиц в кластере 
end
cluster = CL;
% найдем номера частиц вне кластера(дающие вклад в регулярную силу)
for z=1:s
    rcl = setdiff(num_b,CL(z,:));
    k = size(find(rcl),2);
    if k~=0
        for i=1:k
            rCL(z,i)=rcl(i);
        end
    else
        rCL(z) = 0;
    end
%     rCL(z,:) = rcl;
end
out_claster=rCL;
%составим вектора номеров частиц с соседями и без 
bwn=0;% - body with neighbour
bwon=0;% - body without neighbors
for z = 1:s
    bwn =union(CL(z,:),bwn);
end
bwon = setdiff(num_b,bwn);
bwn = setdiff(num_b,bwon);
end