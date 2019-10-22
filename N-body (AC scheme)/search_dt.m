%2-�� ���: ���������� min dt ��� ������ �������:
function [r, t]  = search_dt(T,Nb,M,J,X,Y,Z)

kmax = size(J,2);% - max ����� ������� � �������
alf = 0.5/T;

for i = 1:Nb
    for j = 1:kmax
        R(i,j) = 0;
    end
end
clear kmax;

for i = 1:Nb
    %� ������ ���� � ������� ��� �������
    if J(i,1) == 0
        minR(i) = 0;% - ���������� ����� �������� � � ��������� �������
        dt(i) = 0;% ����������� ���������� ��� dT(n,i)
    %� ������ ���� � ������� ���� ������
    else
        k = size(find(J(i,:)),2);% - max ����� ������� ��� i-�� �������
        minR(i) = 0;% - ���������� ����� �������� � � ��������� �������
        %���������� ���������� ����� i-�� �������� � � ��������
        for j = 1:k 
            R(i,j) = ((X(J(i,j))-X(i))^2+(Y(J(i,j))-Y(i))^2+(Z(J(i,j))-Z(i))^2)^(1/2);   
            if R(i,j)<minR(i) || minR(i)==0
                minR(i)= R(i,j);
            end
        end
        dt(i) = alf*minR(i)^(3/2);%alf*minR(i)^(3/2);%����� min ���������� ��� dt ��� ������ �������
        if minR(i)<0.05 && M(i)<1e1
            dt(i) = dt(i)/3;
        end
    end
end
clear reg; 
r = R;
t = dt;
end