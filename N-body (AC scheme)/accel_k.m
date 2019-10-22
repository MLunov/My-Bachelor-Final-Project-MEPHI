function [Ax, Ay, Az] = accel_k(Nb,M,CL,CL_out,X,Y,Z)
A = 1.393825555; % G*M_sun/[1��.�����.]^2 (��/c^2)
for i = 1:Nb
    ax = zeros(Nb,Nb);
    ay = zeros(Nb,Nb);
    az = zeros(Nb,Nb);
    sx=0; sy=0; sz=0;
    %������� �������������� i-�� ������� � ��������
    if find(CL==i) ~= 0
        [z,d] = find(CL==i);%(1-� ����� - ����� ��������;2-� - ����� ��� � ��������)
        J = size(find(CL_out(z,:)),2);%-����� ��� ��� ��������   
    else
        J = Nb;
    end
    if X(i)==inf && Y(i)==inf && Z(i)==inf
            Ax(i)=0;%��������� i-�� �������
            Ay(i)=0;%� n-�� ������ �������
            Az(i)=0;
    else
        for j = 1:J
            if find(CL==i) ~= 0
                l = CL_out(z,j);%-�������������
            else
                l=j;
            end
            ax(i,j)=0;
            ay(i,j)=0;
            az(i,j)=0;
            if X(l)==inf && Y(l)==inf && Z(l)==inf 
                ax(i,j)=0;
                ay(i,j)=0; 
                az(i,j)=0; 
            else
                if  X(l)-X(i)==0 && Y(l)-Y(i)==0 && Z(l)-Z(i)==0%�������� �� ��, ��� l!=i
                    ax(i,j)=0;
                    ay(i,j)=0; 
                    az(i,j)=0;  
                else
                    %��������, ��� �������������� ����� i � j ���������
                    if ax(j,i)~=0 | ay(j,i)~=0 | az(j,i)~=0 
                        ax(i,j)= -ax(j,i)*(M(l)/M(i)); 
                        ay(i,j)= -ay(j,i)*(M(l)/M(i)); 
                        az(i,j)= -az(j,i)*(M(l)/M(i));
                    else
                        amod=M(l)/((X(l)-X(i))^2+(Y(l)-Y(i))^2+(Z(l)-Z(i))^2);
                        alp=atan(abs((Z(l)-Z(i))/(((Y(l)-Y(i))^2+(X(l)-X(i))^2)^(1/2))));
                        if  X(l)-X(i)==0 && Y(l)-Y(i)==0
                            beta = 0;
                        else
                            beta=atan(abs((Y(l)-Y(i))/(X(l)-X(i))));
                        end
                        xproj=cos(alp)*cos(beta);
                        yproj=cos(alp)*sin(beta);  
                        zproj=sin(alp);
                        ax(i,j)=amod*xproj;
                        ay(i,j)=amod*yproj;
                        az(i,j)=amod*zproj;
                        if X(l)-X(i)>0
                            ax(i,j)=ax(i,j);
                        elseif X(l)-X(i)<0
                            ax(i,j)=-ax(i,j);
                        else ax(i,j)=ax(i,j)*0;
                        end
                        if Y(l)-Y(i)>0
                            ay(i,j)=ay(i,j);
                        elseif Y(l)-Y(i)<0
                            ay(i,j)=-ay(i,j);
                        else ay(i,j)=ay(i,j)*0;
                        end
                        if Z(l)-Z(i)>0
                            az(i,j)=az(i,j);
                        elseif Z(l)-Z(i)<0
                            az(i,j)=-az(i,j);
                        else az(i,j)=az(i,j)*0;
                        end
                    end
                end
            end
            sx=sx+ax(i,j);
            sy=sy+ay(i,j);
            sz=sz+az(i,j);
        end
        Ax(i)=sx*A;%��������� i-�� �������
        Ay(i)=sy*A;%� n-�� ������ ������� ������������� 
        Az(i)=sz*A;%���������� �����
    end
end
end
