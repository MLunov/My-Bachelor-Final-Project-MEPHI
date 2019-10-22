function [SAx,SAy,SAz,X,Y,Z,VX,VY,VZ,M,Rg] = accel_s(bwn,ts,K,s,imin,CL,M,X,Y,Z,VX,VY,VZ,KAx,KAy,KAz)

A = 1.393825555; % G*M_sun/[1ед.расст.]^2 (км/c^2)
rg = 9.56*(10^(-6));% - константа для радиуса Шварцшильда[2.95(км)/1(ед.расст.)]
Nb = size(find(bwn),2);
% Rg = 0.1;
for j = 1:Nb
    i = bwn(j);%-переобозначим
    SX(1,i) = X(i); SY(1,i) = Y(i); SZ(1,i) = Z(i);
    SVX(1,i) = VX(i); SVY(1,i) = VY(i); SVZ(1,i) = VZ(i);
    SAX(1,i) = 0;  SAY(1,i) = 0;  SAZ(1,i) = 0;
    SM(1,i) = M(i);
end
for j = 1:Nb
    p(bwn(j))=0;
end
for z = 1:s
    t = ts(imin(z));% - иррегулярный временной шаг для кластера i
    for k = 1:K(imin(z))
        nb = size(find(CL(z,:)),2);% число тел в кластере i
%         nb = find(CL(z,:));
%         nb = size(nb,2)% число тел в кластере i
        bMAX=0;% - тело с которым происходит наибольшее число столкновений в один шаг по t
        MAX=0;% - маскимальное число столкновений для i-го тела
        ax = zeros(Nb,Nb);
        ay = zeros(Nb,Nb);
        az = zeros(Nb,Nb);
        %обновляем вклад от ИРРЕГУЛЯРНОЙ силы
        for l = 1:nb       
            i = CL(z,l);
%             Rg = rg*SM(k,i)*10000;% - радиус Шварцшильда
            Rg_i = 3*rg*SM(k,i);% - радиус Шварцшильда
            b=0; mas=0; 
            MRx=0; MRy=0; MRz=0;
%             ox=0; oy=0; oz=0; 
            px=0; py=0; pz=0; 
            sx=0; sy=0; sz=0;
            if SX(k,i)==inf && SY(k,i)==inf && SZ(k,i)==inf 
                SAX(k,i)=0;%иррегулярное ускорение i-ой частицы
                SAY(k,i)=0;%в k-ый момент времени
                SAZ(k,i)=0;
                intAx(k,i) = 0;%intermediate
                intAy(k,i) = 0;
                intAz(k,i) = 0;
            else
                for j = 1:nb
                    ax(i,j)=0;
                    ay(i,j)=0;
                    az(i,j)=0; 
                    Rg_j = 3*rg*SM(k,CL(z,j));
                    if CL(z,j) == i
                        ax(i,j)=0;
                        ay(i,j)=0;
                        az(i,j)=0;
                    else
                        if SX(k,CL(z,j))==inf && SY(k,CL(z,j))==inf && SZ(k,CL(z,j))==inf
                            ax(i,j)=0;
                            ay(i,j)=0; 
                            az(i,j)=0;                                                    
                            %условие столкновения тел(З.С.И.)
                        elseif abs(SX(k,CL(z,j))-SX(k,i))<Rg_i+Rg_j && abs(SY(k,CL(z,j))-SY(k,i))<Rg_i+Rg_j && abs(SZ(k,CL(z,j))-SZ(k,i))<Rg_i+Rg_j
                            if bMAX~=0
                                %проверка, что сталкиваются тела одной группы
                                if abs(SX(k,bMAX)-SX(k,i))>=Rgmax+Rg_i || abs(SY(k,bMAX)-SY(k,i))>=Rgmax+Rg_i || abs(SZ(k,bMAX)-SZ(k,i))>=Rgmax+Rg_i
                                    MAX=0;
                                end
                            end
                            px = px + SM(k,CL(z,j))*SVX(k,CL(z,j));
                            py = py + SM(k,CL(z,j))*SVY(k,CL(z,j));
                            pz = pz + SM(k,CL(z,j))*SVZ(k,CL(z,j));
                            MRx = MRx+SM(k,CL(z,j))*(SX(k,CL(z,j)));
                            MRy = MRy+SM(k,CL(z,j))*(SY(k,CL(z,j)));
                            MRz = MRz+SM(k,CL(z,j))*(SZ(k,CL(z,j)));
%                             ox = ox + SX(k,CL(z,j));
%                             oy = oy + SY(k,CL(z,j));
%                             oz = oz + SZ(k,CL(z,j));
                            mas = mas + SM(k,CL(z,j));
                            SX(k+1,CL(z,j))=inf; SVX(k+1,CL(z,j))=0; 
                            SY(k+1,CL(z,j))=inf; SVY(k+1,CL(z,j))=0;
                            SZ(k+1,CL(z,j))=inf; SVZ(k+1,CL(z,j))=0;
                            SM(k+1,CL(z,j))=0;
                            p(CL(z,j))=1;
                            ax(i,j)=0; ay(i,j)=0; az(i,j)=0;
                            b=b+1; % - счетчик столкновений 
                            if b>=MAX
                                MAX=b;
                                bMAX=i; % тело с которым происходит max столкновений
                                rrg(l) = Rg_i;
                                Rgmax = Rg_i;
                                PX=px;PY=py;PZ=pz;
                                MRX = MRx; MRY = MRy; MRZ = MRz;
%                                 OX=ox;OY=oy;OZ=oz;
                                mmax=mas;
                            end
                        else
                            amod=SM(k,CL(z,j))/((SX(k,CL(z,j))-SX(k,i))^2+(SY(k,CL(z,j))-SY(k,i))^2+(SZ(k,CL(z,j))-SZ(k,i))^2);
                            alp=atan(abs((SZ(k,CL(z,j))-SZ(k,i))/(((SY(k,CL(z,j))-SY(k,i))^2+(SX(k,CL(z,j))-SX(k,i))^2)^(1/2))));
                            if  SX(k,CL(z,j))-SX(k,i)==0 && SY(k,CL(z,j))-SY(k,i)==0
                                beta = 0;
                            else
                                beta=atan(abs((SY(k,CL(z,j))-SY(k,i))/(SX(k,CL(z,j))-SX(k,i))));
                            end
                            xproj=cos(alp)*cos(beta);
                            yproj=cos(alp)*sin(beta);        
                            zproj=sin(alp);
                            ax(i,j)=amod*xproj;
                            ay(i,j)=amod*yproj;
                            az(i,j)=amod*zproj;
                            if SX(k,CL(z,j))-SX(k,i)>0
                                ax(i,j)=ax(i,j);
                            elseif SX(k,CL(z,j))-SX(k,i)<0
                                ax(i,j)=-ax(i,j);
                            else ax(i,j)=ax(i,j)*0;
                            end
                            if SY(k,CL(z,j))-SY(k,i)>0
                                ay(i,j)=ay(i,j);
                            elseif SY(k,CL(z,j))-SY(k,i)<0
                                ay(i,j)=-ay(i,j);
                            else ay(i,j)=ay(i,j)*0;
                            end
                            if SZ(k,CL(z,j))-SZ(k,i)>0
                                az(i,j)=az(i,j);
                            elseif SZ(k,CL(z,j))-SZ(k,i)<0
                                az(i,j)=-az(i,j);
                            else az(i,j)=az(i,j)*0;
                            end
                        end
                    end
                    sx=sx+ax(i,j);
                    sy=sy+ay(i,j);
                    sz=sz+az(i,j);
                end
                SAX(k,i)=sx*A;%ускорение i-ой частицы
                SAY(k,i)=sy*A;%в k-ый момент времени приобретенное 
                SAZ(k,i)=sz*A;%иррегулярной силой(от соседних частиц)
                % вычислим суммарное ускорение приобретенное i-ой частицей
                intAx(k,i) = KAx(i) + SAX(k,i);%intermediate
                intAy(k,i) = KAy(i) + SAY(k,i);
                intAz(k,i) = KAz(i) + SAZ(k,i);
                if k == 1
%                     if SAX(k,i)==0 && SAY(k,i)==0 && SAZ(k,i)==0 
%                         SAx(i)=0;
%                         SAy(i)=0;
%                         SAz(i)=0;
%                     else
                    SAx(i)=SAX(k,i);%intermediate
                    SAy(i)=SAY(k,i);
                    SAz(i)=SAZ(k,i);
%                     Ax(i) = intAx(k,i);
%                     Ay(i) = intAy(k,i);
%                     Az(i) = intAy(k,i); 
                end
            end
            % Вычисление координат и скоростей тел через t
            if b>0 % в случае столкновения высчитываем новые r и V через З.С.И.
                sumM=mmax+SM(k,bMAX);
                MRxmax = SM(k,bMAX)*(SX(k,bMAX));
                MRymax = SM(k,bMAX)*(SY(k,bMAX));
                MRzmax = SM(k,bMAX)*(SZ(k,bMAX));
                sumMRX = MRX+MRxmax;
                sumMRY = MRY+MRymax;
                sumMRZ = MRZ+MRzmax;
                SX(k+1,bMAX) = sumMRX/sumM;
                SY(k+1,bMAX) = sumMRY/sumM;
                SZ(k+1,bMAX) = sumMRZ/sumM;
%                 Xmean = (OX+SX(k,bMAX))/(MAX+1);
%                 Ymean = (OY+SY(k,bMAX))/(MAX+1);
%                 Zmean = (OZ+SZ(k,bMAX))/(MAX+1);
%                 if SM(k,bMAX)/mmax == 1
%                     SX(k+1,bMAX) = Xmean;
%                     SY(k+1,bMAX) = Ymean;
%                     SZ(k+1,bMAX) = Zmean;
%                 elseif SM(k,bMAX)/mmax > 1
%                     KF = koeff(SM(k,bMAX),mmax);
%                     SX(k+1,bMAX) = Xmean + (SX(k,bMAX) - Xmean)*KF;
%                     SY(k+1,bMAX) = Ymean + (SY(k,bMAX) - Ymean)*KF;
%                     SZ(k+1,bMAX) = Zmean + (SZ(k,bMAX) - Zmean)*KF;
%                 else
%                     KF = koeff(mmax,SM(k,bMAX));
%                     SX(k+1,bMAX) = Xmean + (OX/MAX - Xmean)*KF;
%                     SY(k+1,bMAX) = Ymean + (OY/MAX - Ymean)*KF;
%                     SZ(k+1,bMAX) = Zmean + (OZ/MAX - Zmean)*KF;                    
%                 end
                SVX(k+1,bMAX)= (PX+SM(k,bMAX)*SVX(k,bMAX))/sumM;
                SVY(k+1,bMAX)= (PY+SM(k,bMAX)*SVY(k,bMAX))/sumM;
                SVZ(k+1,bMAX)= (PZ+SM(k,bMAX)*SVZ(k,bMAX))/sumM;
                SM(k+1,bMAX)=SM(k,bMAX)+mmax;          
            elseif p(i)~=0 && Rg_i<max(rrg)
                SX(k+1,i)=inf; SVX(k+1,i)=0; 
                SY(k+1,i)=inf; SVY(k+1,i)=0;
                SZ(k+1,i)=inf; SVZ(k+1,i)=0;
                SM(k+1,i)=0;
            else
                SX(k+1,i)=SX(k,i)+(t*SVX(k,i)+(t^2)*intAx(k,i)/2)/(3.08568*10^5);
                SY(k+1,i)=SY(k,i)+(t*SVY(k,i)+(t^2)*intAy(k,i)/2)/(3.08568*10^5);
                SZ(k+1,i)=SZ(k,i)+(t*SVZ(k,i)+(t^2)*intAz(k,i)/2)/(3.08568*10^5);
                SVX(k+1,i)=SVX(k,i)+t*intAx(k,i);
                SVY(k+1,i)=SVY(k,i)+t*intAy(k,i);
                SVZ(k+1,i)=SVZ(k,i)+t*intAz(k,i);
                SM(k+1,i)=SM(k,i);  
            end
        end
    end
    %Получим координаты и скорости тел через T
    for l = 1:nb        
        i = CL(z,l); 
        X(i)=SX(k+1,i);
        Y(i)=SY(k+1,i);
        Z(i)=SZ(k+1,i);
        VX(i)=SVX(k+1,i);
        VY(i)=SVY(k+1,i);
        VZ(i)=SVZ(k+1,i);
        M(i)=SM(k+1,i);
    end
end
for i=1:Nb
    Rg(i) = 3*rg*SM(1,bwn(i));
end
end