clc
close all; 
clear all; close all;
%% Definizione profilo
% naca = 'NACA0002';
% n = 100;
% [xp,yp] = NACA_generator(naca,n,'cos',0); grid on;

% naca = load('64212.mat');
% dummy = size(naca.air,1); dummy = [1,ceil(linspace(2,dummy-1,n*2-3)),dummy];
% xp = naca.air(dummy,1)';
% yp = naca.air(dummy,2)';

dummy = load('acu.mat');
xp = dummy.xp2;
yp = dummy.yp2;
dummy = load('antiorario.dat');
xp = dummy(:,1)';
yp = dummy(:,2)';

pan = [xp;yp];

Npt = length(xp);
Npan = length(xp)-1; % numero pannelli

alpha_v = [-2.5:0.5:2.5]%:20]%:2:18];

for av = 1:size(alpha_v,2)
alpha_v(av)

alpha = alpha_v(av)*pi/180; % incidenza profilo
U_inf = 135*[cos(alpha) sin(alpha)]; % velocità asintotica
%% Punti di collocazione
% [lfun] = HS_staz(xp,yp,U_inf);
% collocation_points = lfun{1};

%[xc,yc] = collocation_points(xp,yp);



%title([naca(1:4) ' ' naca(5:end)])
%hold off
[~,A,B,xc,yc,dl,N,T,x] = HS_staz(xp,yp,U_inf);
U_inf_t = U_inf*T;
u_pert = B*x + U_inf_t';



% figure
% plot(xc,x(1:end-1))
%% Coefficiente di pressione
cp = 1-(u_pert./norm(U_inf)).^2;
DCp(av) = min(cp) - cp(1);

figure(av*3 -2)
axis equal
hold on
plot(xp,yp,'Color',[0 0.5 1],'LineWidth',2);
plot(xc,yc,'xr','LineWidth',2);


figure(3*av-1)
n = length(xc)/2; 
hold on
plot(xc(1:n),cp(1:n),'Color',[0 0.5 1],'LineWidth',2);
plot(xc(n:end),cp(n:end),'Color',[1 0.5 0],'LineWidth',2);
grid on
set(gca,'YDir','Reverse'); % reverse y axis
alfa=num2str(alpha*180/pi);
tit1 = 'COEFFICIENTE DI PRESSIONE';
tit2 = ['Incidenza: ',alfa, ' gradi'];
title({tit1;tit2})
legend('dwn','up')
xlabel('corda');
ylabel('-C_p');




%***********************************************************************************
%% confronto con xfoil automatizzato
%prof=str2double(naca(5:end));
%formatSpec = 'naca %2d \n oper.i \n mach %1.5f \n alfa %2.1f \n cpwr cp%2d.dat \n';
%id=fopen('xf.txt','w');
%fprintf(id,formatSpec,prof,norm(U_inf)/343,alpha*180/pi,prof);
%fclose(id);

%!D:\Documenti\Polimi\xfoil\xfoil.exe<xf.txt
%prof=num2str(prof);
%[x_xfoil,y_xfoil,cp_xfoil] = xfoilimport(['cp' prof '.dat']);

%figure(3)
%plot(x_xfoil,cp_xfoil,'g*')
end

%plot(alpha_v,abs(DCp))