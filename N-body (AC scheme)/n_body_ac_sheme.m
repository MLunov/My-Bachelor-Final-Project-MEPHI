clc
clear all

%1 ��. ���������  = 3.08568*10^5 �� = 1*10^(-8) �� 
%1 ��.����� = M_sun; m = M_pbh/M_sun
%1 ��.������� = 1 ���.
% G*M_sun = 132712*10^6;%(��^3/c^2)
% G*M_sun = 132712438*10^12; %(�^3/c^2)
% A = 1,393825555;% G*M_sun/[1��.�����.]^2 (��/c^2) 
% B = 430089.963962562546991;% G*M_sun/[1��.�����.] (��/c)^2 
% B2 = 4.3009e11;% G*M_sun/[1��.�����.] (�/c)^2 
% rg= 9,56*10^(-6);% - ��������� ��� ������� �����������[2.95(��)/1(��.�����.)]

%--------------------------------------------------------------------------
x = [5e2]; y = [5e2]; z = [5e2];
Nb1 = length(x);%����� ���
for i=1:Nb1 %�������� 1-�� ���� �� t ������� ��������� ���
    X(1,i) = x(i); Y(1,i) = y(i); Z(1,i) = z(i);
end
Vx = [0]; Vy = [0]; Vz = [0];
for i=1:Nb1 %�������� 1-�� ���� �� t ������� ��������� ���
    VX(1,i) = Vx(i); VY(1,i) = Vy(i); VZ(1,i) = Vz(i);
end
m = [1e4];
for i=1:Nb1 %�������� 1-�� ���� �� t ������� ���� ���
    M(1,i) = m(i);
end

% Nb2=5;Nb3=13;Nb4=29;Nb=50;
Nb=50;

% mmax = 1e4; m = mmax*rand(1,Nb2);
% for i=(Nb1+1):Nb2 %�������� 1-�� ���� �� t ������� ���� ���
%     M(1,i) = m(i);
% end
% 
% mmax = 1e3; m = mmax*rand(1,Nb3);
% for i=(Nb2+1):Nb3 %�������� 1-�� ���� �� t ������� ���� ���
%     M(1,i) = m(i);
% end
% 
% mmax = 1e2; m = mmax*rand(1,Nb4);
% for i=(Nb3+1):Nb4 %�������� 1-�� ���� �� t ������� ���� ���
%     M(1,i) = m(i);
% end

mmax = 1e1; m = mmax*rand(1,Nb);
for i=(Nb1+1):Nb %�������� 1-�� ���� �� t ������� ���� ���
    M(1,i) = m(i);
end

sumM=0;
for i = 1:Nb
    sumM = sumM + M(1,i);
end

ROS = 1e3;% radius of system
Xmax = ROS; x = Xmax*rand(1,Nb);
Ymax = ROS; y = Ymax*rand(1,Nb);
Zmax = ROS; z = Zmax*rand(1,Nb);
clear Xmax Ymax Zmax;
for i=(Nb1+1):Nb %�������� 1-�� ���� �� t ������� ��������� ���
    X(1,i) = x(i); Y(1,i) = y(i); Z(1,i) = z(i);
end

B = 430089.963962562546991;% G*M_sun/[1��.�����.] (��/c) 
% V_virial = sqrt((sumM*B)/ROS);% ������� �������� �������
% clear B;

% Vx = -V_virial+2*V_virial*rand(1,Nb);
% Vy = -V_virial+2*V_virial*rand(1,Nb);
% Vz = -V_virial+2*V_virial*rand(1,Nb);
Vx = 0*ones(1,Nb);Vy = 0*ones(1,Nb);Vz = 0*ones(1,Nb);
for i=2:Nb %�������� 1-�� ���� �� t ������� ��������� ���
    VX(1,i) = Vx(i); VY(1,i) = Vy(i); VZ(1,i) = Vz(i);
end

%������ ������� ���:
for i = 1:Nb
    num_b(i)=i;
end
clear mmax;
%--------------------------------------------------------------------------

%������� ���������� ��� � 3-� ������������ � ������ t = 0
%1 ��. ���������  = 3.08568*10^5 �� = 1*10^(-8) �� (���������� ����������
%������ ���� ������ ��� �� 1 �������)

% Nb = 50;
% Xmax = 20; x = Xmax*rand(1,Nb);
% Ymax = 20; y = Ymax*rand(1,Nb);
% Zmax = 20; z = Zmax*rand(1,Nb);
% clear Xmax; clear Ymax; clear Zmax;

% x = [0 100];
% y = [0 0];
% z = [100 100];

% Nb = length(x);%����� ���

%������� �� �������� � 3-� ������������ � ������ t = 0
%1 ��. ��������  = 1 ��/�
% Vx = 0*ones(1,Nb);Vy = 0*ones(1,Nb);Vz = 0*ones(1,Nb);

% ������� �� ����� (1��.����� = M_sun; m = M_pbh/M_sun)
% m = 1*ones(1,Nb);
% m = [1e3 1e3];


% for i=1:Nb %�������� 1-�� ���� �� t ������� ��������� ���
%     X(1,i) = x(i); Y(1,i) = y(i); Z(1,i) = z(i);
% end
% 
% for i=1:Nb %�������� 1-�� ���� �� t ������� ��������� ���
%     VX(1,i) = Vx(i); VY(1,i) = Vy(i); VZ(1,i) = Vz(i);
% end
% 
% for i=1:Nb %�������� 1-�� ���� �� t ������� ���� ���
%     M(1,i) = m(i);
% end
% 
% %������ ������� ���:
%     for i = 1:Nb
%         num_b(i)=i;
%     end
%     
clear x; clear y; clear z; 
clear Vx; clear Vy; clear Vz;
clear m; 
clear i;

%������� ������� �������
t_0=0;
t_n=2e4;
T = 1; % - ������ ���������� ��� ������� (1 ��.������� = 1 ���.)
hT = t_0:T:t_n; % - ������ �������� ����� �� �������
N = length(hT); % - ������ ����� ��������� �����(N)
clear t_0 t_n hT;

% N=1;
for n=1:N
    %1-�� ���: ����� ������� ��� ������ �������
    [Rn,J] = select(Nb,M(n,:),X(n,:),Y(n,:),Z(n,:));
    Ros(n,:)=Rn;% - ������ ������ �������(Radius of selection)
    J;% - ������� ������� ������� ��� ������ �������
    clear Rn;
    
    %2-�� ���: ��������, ���� �� �����(�) ���� � ����� �������:
    if find(J)~=0
        %���� ��:
        %��������� ���� �� ���������� � ������������ ������������
        
        %3-�� ���: ���������� dt ��� ������ ������� � ������� ���� �����:
        [R, t] = search_dt(T,Nb,M(n,:),J,X(n,:),Y(n,:),Z(n,:));
%         R;% - ������� ���������� ����� �������� � � ��������
        dt(n,:)=t;% - ������ min ���������� ����� dt ��� ������ ������� �� ������ T
        clear t;
        
        %4-�� ���: ����� ������ � min dt(������������ ��������� �����)
        %� ����������� �� ��������� �������
        [im, dtm, s, CL,CL_out,bw,bwo] = mindt(dt(n,:),M(n,:),J,num_b);
%         s;% - ����� ����� ������
        for i = 1:s
            imin(n,i)=im(i);% - ������ ����� ������
            dtmin(n,i)=dtm(i);% - ������������ ��������� ��� ��� imin ������ � �� ���������
        end
%         Cl;% - ������� ������� ������� � �������� imin �������
%         CL;% - ������� ������� ������ � imin-��������(� ������ imin �������)
%         CL_out;% - ������� ������� ������ ��� imin-�� ��������(������ ����� � ���������� ����)
        k = size(find(bw),2);
        for i = 1:k
            bwn(n,i) = bw(i);% - ������ ������� ���� ������ � ��������(body with neighbour)
        end
        if k~=Nb
            for i = 1:(Nb-k)
                bwon(n,i) = bwo(i);% - ������ ������� ���� ������ ��� �������(body without neighbors)
            end
        else
            bwon(n,:)=0;
        end
        clear im dtm bw bwo;
        
        %5-�� ���: ���������� R � V ���� ������ ������ dT
        %(��� ��������� imin �� ��� ��� ���� (����� dt) = dT)
        %������� ���������� ��������� ��� ���� ������ 
        if find(CL_out) ~= 0
            [kax, kay, kaz] = accel_k(Nb,M(n,:),CL,CL_out,X(n,:),Y(n,:),Z(n,:));            
            KAx(n,:)=kax; KAy(n,:)=kay; KAz(n,:)=kaz;
            clear kax, clear kay, clear kaz;
            %������� R � V ��� ��� ��������� 
            k = size(find(bwon(n,:)),2);
            for j = 1:k
                i = bwon(n,j);%-�������������
                Ax(n,i) = KAx(n,i);
                Ay(n,i) = KAy(n,i);
                Az(n,i) = KAz(n,i);
                %��������� ���������� � �������� ��� ����� T
                X(n+1,i)=X(n,i)+(T*VX(n,i)+(T^2)*Ax(n,i)/2)/(3.08568*10^5);
                Y(n+1,i)=Y(n,i)+(T*VY(n,i)+(T^2)*Ay(n,i)/2)/(3.08568*10^5);
                Z(n+1,i)=Z(n,i)+(T*VZ(n,i)+(T^2)*Az(n,i)/2)/(3.08568*10^5);
                VX(n+1,i)=VX(n,i)+T*Ax(n,i);%km
                VY(n+1,i)=VY(n,i)+T*Ay(n,i);
                VZ(n+1,i)=VZ(n,i)+T*Az(n,i);
                M(n+1,i)=M(n,i);
            end
        else
            for j = 1:Nb
                KAx(n,i)=0;
                KAy(n,i)=0;
                KAz(n,i)=0;
            end  
        end
        clear i j k;
        %������������� �������: ������� ����� ������� ����� dt ��� ������
        %������� � mindt � � ��������
        for z = 1:s
            K(n,imin(n,z)) = ceil(T/dtmin(n,z));% - ����� ������������ ����� ��� ������� � dtmin
            dts(n,imin(n,z)) = T/(ceil(T/dtmin(n,z)));% - �� dtmin ��� ����������
        end
        clear z;
        %������� ������������ ��������� ��� ������ � ���������
        %� ������ ������� �� R � V  
        [sax,say,saz,SX,SY,SZ,SVX,SVY,SVZ,SM,rg] = accel_s(bwn(n,:),dts(n,:),K(n,:),s,imin(n,:),CL,M(n,:),X(n,:),Y(n,:),Z(n,:),VX(n,:),VY(n,:),VZ(n,:),KAx(n,:),KAy(n,:),KAz(n,:));
        for i = 1:Nb
            if find(bwn(n,:)==i)~=0
                SAx(n,i)=sax(i);
                SAy(n,i)=say(i);
                SAz(n,i)=saz(i);
            else
                SAx(n,i)=0;
                SAy(n,i)=0;
                SAz(n,i)=0;
            end
        end
        %������� R � V ��� � ���������
        k = size(find(bwn(n,:)),2);
        for j = 1:k   
            i = bwn(n,j);%-�������������
            Ax(n,i) = KAx(n,i)+SAx(n,i);
            Ay(n,i) = KAy(n,i)+SAy(n,i);
            Az(n,i) = KAz(n,i)+SAz(n,i);
            %���������� � �������� ��� ����� T
            X(n+1,i)=SX(i);
            Y(n+1,i)=SY(i);
            Z(n+1,i)=SZ(i);
            VX(n+1,i)=SVX(i);
            VY(n+1,i)=SVY(i);
            VZ(n+1,i)=SVZ(i);
            M(n+1,i)=SM(i);
        end
        for i = 1:k
            Rg(n,i) = rg(i);% -�������������� ������ ���
        end
        
        clear sax say saz; clear SX SY SZ; clear SVX SVY SVZ; clear SM;
        clear rg; clear i j k; 
    else %���� ���:
        %���������� ������ ���������� ������������
        [ax, ay, az] = accel(Nb,M(n,:),X(n,:),Y(n,:),Z(n,:)); 
        Ax(n,:)=ax; Ay(n,:)=ay; Az(n,:)=az;
        clear ax, clear ay, clear az;
        %�������� R � V ��� 
        for i = 1:Nb
            %��������� ���������� � �������� ��� ����� T
            X(n+1,i)=X(n,i)+(T*VX(n,i)+(T^2)*Ax(n,i)/2)/(3.08568*10^5);
            Y(n+1,i)=Y(n,i)+(T*VY(n,i)+(T^2)*Ay(n,i)/2)/(3.08568*10^5);
            Z(n+1,i)=Z(n,i)+(T*VZ(n,i)+(T^2)*Az(n,i)/2)/(3.08568*10^5);
            VX(n+1,i)=VX(n,i)+T*Ax(n,i);%km
            VY(n+1,i)=VY(n,i)+T*Ay(n,i);
            VZ(n+1,i)=VZ(n,i)+T*Az(n,i);
            M(n+1,i)=M(n,i);          
        end
    end
end
for n =1:N
    Nb_n(n) = size(find(M(n,:)),2);
end
clear hh NN dNb
hh = 1000
NN = 1:hh:N
for n = 1:(length(NN)-1)
    dNb(n) = abs((Nb_n(NN(n+1))-Nb_n(NN(n)))/(hh*T));
end
% % ��� �������� ������ �������� ������� ���� � ������ ��� �������
% for n=1:N
%     for i = 1:Nb
%         Vmod(n,i) = sqrt((VX(n,i))^2+(VY(n,i))^2+(VZ(n,i))^2);
%     end
% end
% % ���������� ������������ ������� ������� � n-�� ������ �������
% for n = 1:N
%     Et = 0;
%     for i = 1:Nb
%         Vmod(n,i) = sqrt((VX(n,i))^2+(VY(n,i))^2+(VZ(n,i))^2);
%         Et = Et+M(n,i)*(Vmod(n,i))^2;  
%     end
%     E_t(n) = (Et/2)/1e6;
% end
% % ���������� ������������� ������� ������� � 1-�� ������ �������
% B = 430089.963962562546991;% G*M_sun/[1��.�����.] (��/c)^2
% n = 1; Eu=0;
% for i =1:(Nb-1)
%     for j = (i+1):Nb
%         Eu = Eu+(M(n,j)*M(n,i))/(sqrt((X(n,j)-X(n,i))^2+(Y(n,j)-Y(n,i))^2+(Z(n,j)-Z(n,i))^2));
%     end
%     U_0 = B*Eu;
% end
% % ���������� ������������� ������� ������� � n-�� ������ �������
% for n = 1:N
%     Eu=0;
%     for i =1:(Nb-1)
%         for j = (i+1):Nb
%             Eu = Eu+(M(n,j)*M(n,i))/(sqrt((X(n,j)-X(n,i))^2+(Y(n,j)-Y(n,i))^2+(Z(n,j)-Z(n,i))^2));
%         end
%     end
%     if Eu == 0
%         break
%     end
%     E_u(n) = (U_0-B*Eu)/1e6;
%     E_total(n) = E_t(n)+E_u(n);
% end
% clear Et Eu i j;

%%
% v = VideoWriter('check.avi');
% open(v);
d=10;
h=50;
% Nb1=1;Nb2=5;Nb3=12;Nb4=25;Nb5=50;Nb=100;
% Nb1=1;Nb2=3;Nb3=7;Nb4=13;Nb5=25;Nb=50;
% Nb1 =1;Nb2=5;Nb3=13;Nb4=29;Nb=50;
for n=1:h:N
    subplot(1,2,1);
    plot3(X(n,:),Y(n,:),Z(n,:),'ko')
    hold on
%     plot3(X(n,1),Y(n,1),Z(n,1),'r*')
%     for i = (Nb1+1):Nb2
%     plot3(X(n,i),Y(n,i),Z(n,i),'y*')
%     end
%     for i = (Nb2+1):Nb3
%     plot3(X(n,i),Y(n,i),Z(n,i),'g*')
%     end
%     for i = (Nb3+1):Nb4
%     plot3(X(n,i),Y(n,i),Z(n,i),'c*')
%     end
%     for i = (Nb4+1):Nb5
%     plot3(X(n,i),Y(n,i),Z(n,i),'b*')
%     end
%     for i = (Nb5+1):Nb
%     plot3(X(n,i),Y(n,i),Z(n,i),'k*')
%     end
% %--------------------------------------------------------------------------
% %     plot3(X(n,16),Y(n,16),Z(n,16),'r+')
% %     legend('~{10}^6','~{10}^5','~{10}^4','~{10}^3','~{10}^2','~{10}^1')
% %     plot(X(n,50),Y(n,50),'o','color','r')
% %--------------------------------------------------------------------------
    grid on
    title('Gravity dynamics of PBH','fontsize',16,'fontname','times','interpreter','latex')
%     axis([min(X(1,:))-d max(X(1,:))+d min(Y(1,:))-d max(Y(1,:))+d min(Z(1,:))-d max(Z(1,:))+d])
%     axis([-5 15 -5 15 -5 15])   
    axis([0 ROS+d 0 ROS+d 0 ROS+d])
    xlabel('1 u.m. = $10^{-8}$ pc','fontsize',14,'interpreter','latex');
    ylabel('1 u.m. = $10^{-8}$ pc','fontsize',14,'interpreter','latex');
    zlabel('1 u.m. = $10^{-8}$ pc','fontsize',14,'interpreter','latex'); 
%     text(min(X(1,:))-d,max(Y(1,:))+d,min(Z(1,:))-d,strcat('Time = ',num2str((n-1)*T)));
%     text(min(X(1,:))-d,max(Y(1,:))+d,min(Z(1,:))-d-2,strcat('T = ',num2str(T)));
    text(-500,ROS-d,-d,strcat('Time = ',num2str((n-1)*T), ' s'), 'FontSize',14,'fontname','times','interpreter','latex');
    text(-500,ROS-d,-d-50,strcat('$\Delta T$ = ',num2str(T), ' s'), 'FontSize',14,'fontname','times','interpreter','latex');
%     text(-310,ROS-d,-50,'(second)');
    hold off
    getframe;
    
    subplot(1,2,2);
    ii = find(M(n,:)~=0);
    hist(log10(M(n,ii)))
%     hist(M(n,ii),max(M(n,:)/1000))
    grid on
    title('The distribution of $N_b$ from masses','fontsize',16,'fontname','times','interpreter','latex');
%     xlabel('log[M(pbh)/M(sun)]');
    xlabel('$\lg(M_{\rm PBH}/M_\odot)$','fontsize',16,'interpreter','latex');
    ylabel('$N_b$','fontsize',16,'interpreter','latex');
    set(0,'DefaultAxesFontSize',10,'DefaultAxesFontName','Arial Cyr'); 
    getframe;

%     subplot(2,2,4);
%     plot(1:N,Nb_n(:),'linewidth',2)
%     hold on;
%     grid on;
%     title('Dependence of the number of PBH on time','fontsize',16,'fontname','times','interpreter','latex');
%     set(gca,'GridLineStyle','-.','MinorGridLineStyle',':','Layer','top');
%     set(gca,'fontsize',16,'FontName','times');
%     xlabel('$t$','fontsize',16,'interpreter','latex');
%     ylabel('$\Sigma N_i$','fontsize',16,'interpreter','latex');
% %     ylabel('$\frac{N}{\dot{N}}\Delta t$','fontsize',16,'interpreter','latex');
%     xlim([0, 4e4]);
%     %������ �����:
%     set(gcf, 'Position', get(0,'ScreenSize'));
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end
% close(v);
clear Nb1 Nb2 Nb3 Nb4 Nb5;
clear n; clear d; clear i; clear ii;

%��� �������� ���������� ����� �������
% for n=1:N
%     Ms=0;
%     for i = 1:Nb
%         Ms = Ms + M(n,i);
%     end
%     SUMM(n)=Ms;
% end
% clear Ms;

% % ��� �������� ������ �������� ������� ���� � ������ ��� �������
% for n=1:N
%     for i = 1:Nb
%         Vmod(n,i) = sqrt(VX(n,i)^2+VY(n,i)^2+VZ(n,i)^2);
%     end
% end
% 
% % ���������� ������������ ������� ������� � n-�� ������ �������
% for n = 1:N
%     Et = 0;
%     for i = 1:Nb
%         Et = Et+M(n,i)*(Vmod(n,i))^2;
%     end
%     E_t(n) = Et/2;
% end

% % ���������� ������������� ������� ������� � n-�� ������ �������
% for n = 1:N
%     Eu=0;
%     for i =1:(Nb-1)
%         for j = (i+1):Nb
%             Eu = Eu+(M(n,j)*M(n,i))/(sqrt((X(n,j)-X(n,i))^2+(Y(n,j)-Y(n,i))^2+(Z(n,j)-Z(n,i))^2));
%         end
%     end
%     E_u(n) = B*Eu/2;
%     E_total(n) = E_t(n)-E_u(n);
% end

% B = 430089.963962562546991;% G*M_sun/[1��.�����.] (��/c)^2
% for n = 1:N
%     Et = 0;Eu=0;
%     for i =1:(Nb-1)
%         Vmod(n,i) = sqrt(VX(n,i)^2+VY(n,i)^2+VZ(n,i)^2);
%         Et = Et+M(n,i)*(Vmod(n,i))^2;
%         for j = (i+1):Nb
%             Eu = Eu+(M(n,j)*M(n,i))/(sqrt((X(n,j)-X(n,i))^2+(Y(n,j)-Y(n,i))^2+(Z(n,j)-Z(n,i))^2));
%         end
%     end
%     E_t(n) = Et/2;
%     E_u(n) = B*Eu/2;
%     E_total(n) = E_t(n)-E_u(n);
% end
% clear Et Eu i j;

%%
% ��������� ����� ������ �� ������� 
% clear hh NN dNb TT
% hh = 1500
% NN = 1:hh:N
% for n = 1:(length(NN)-1)
%     dNb(n) = abs((Nb_n(NN(n+1))-Nb_n(NN(n)))/(hh*T));
%     TT(n) = Nb_n(NN(n))./dNb(n);
% end
% %
% TT((find(TT==inf))) = 1e5;
% plot(NN(1:end-1),TT,'linewidth',2)
% grid on
% ylim([0,1e4])
% xlim([0,1e4])
% set(gca,'fontsize',24,'FontName','times');
% xlabel('$t$','fontsize',24,'interpreter','latex');
% ylabel('$\frac{N_i}{\dot{N_i}}$','fontsize',28,'interpreter','latex');