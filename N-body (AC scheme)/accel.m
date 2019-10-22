function [Ax, Ay, Az] = accel(Nb,M,X,Y,Z)
A = 1.393825555; % G*M_sun/[1ед.расст.]^2 (км/c^2)
for i = 1:Nb
    ax = zeros(Nb,Nb);
    ay = zeros(Nb,Nb);
    az = zeros(Nb,Nb);
    sx=0; sy=0; sz=0;
    if X(i)==inf && Y(i)==inf && Z(i)==inf
            Ax(i)=0;%ускорение i-ой частицы
            Ay(i)=0;%в n-ый момент времени
            Az(i)=0;
    else
        for j = 1:Nb
            ax(i,j)=0;
            ay(i,j)=0;
            az(i,j)=0;
            if X(j)==inf && Y(j)==inf && Z(j)==inf 
                ax(i,j)=0;
                ay(i,j)=0; 
                az(i,j)=0; 
            else
                if  X(j)-X(i)==0 && Y(j)-Y(i)==0 && Z(j)-Z(i)==0%проверка на то, что j!=i
                    ax(i,j)=0;
                    ay(i,j)=0; 
                    az(i,j)=0;  
                else
                    %проверка, что взаимодействие между i и j посчитано
                    if ax(j,i)~=0 | ay(j,i)~=0 | az(j,i)~=0 
                        ax(i,j)= -ax(j,i)*(M(j)/M(i)); 
                        ay(i,j)= -ay(j,i)*(M(j)/M(i)); 
                        az(i,j)= -az(j,i)*(M(j)/M(i));
                    else
                        amod=M(j)/((X(j)-X(i))^2+(Y(j)-Y(i))^2+(Z(j)-Z(i))^2);
                        alp=atan(abs((Z(j)-Z(i))/(((Y(j)-Y(i))^2+(X(j)-X(i))^2)^(1/2))));
                        if  X(j)-X(i)==0 && Y(j)-Y(i)==0
                            beta = 0;
                        else
                            beta=atan(abs((Y(j)-Y(i))/(X(j)-X(i))));
                        end
                        xproj=cos(alp)*cos(beta);
                        yproj=cos(alp)*sin(beta);  
                        zproj=sin(alp);
                        ax(i,j)=amod*xproj;
                        ay(i,j)=amod*yproj;
                        az(i,j)=amod*zproj;
                        if X(j)-X(i)>0
                            ax(i,j)=ax(i,j);
                        elseif X(j)-X(i)<0
                            ax(i,j)=-ax(i,j);
                        else ax(i,j)=ax(i,j)*0;
                        end
                        if Y(j)-Y(i)>0
                            ay(i,j)=ay(i,j);
                        elseif Y(j)-Y(i)<0
                            ay(i,j)=-ay(i,j);
                        else ay(i,j)=ay(i,j)*0;
                        end
                        if Z(j)-Z(i)>0
                            az(i,j)=az(i,j);
                        elseif Z(j)-Z(i)<0
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
        Ax(i)=sx*A;%ускорение i-ой частицы
        Ay(i)=sy*A;%в n-ый момент времени приобретенное 
        Az(i)=sz*A;%регулярной силой
    end
end
end
