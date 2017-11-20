clc
clear all;
close all;

addpath('./Routines');
%% Definizione profilo
naca = 'NACA0012';
n = 160;
[xp,yp] = NACA_generator(naca,n,'cos',0);

% naca = load('64212.mat');
% dummy = size(naca.air,1); dummy = [1,ceil(linspace(2,dummy-1,n*2-3)),dummy];
% xp = naca.air(dummy,1)';
% yp = naca.air(dummy,2)';

native = load('./nativeNACA0012.txt');

% fig1 = figure(1);
% plot(xp,yp,'bo-');
% title('Airfoil'); grid on; hold on;
% plot(native(:,1),native(:,2),'go-');
% legend('NACA gen','NACA xfoil');

pan = [xp;yp];

Npt  = length(xp);
Npan = length(xp)-1; % numero pannelli

U_mag   = 135;

iter_number = 5000;

SBOiter_max = 1;

% Start with the default options
options = optimoptions('fmincon');
% Modify options setting
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxFunctionEvaluations', 50);
options = optimoptions(options,'OptimalityTolerance', 0.005);
options = optimoptions(options,'FunctionTolerance',   0.005);
options = optimoptions(options,'StepTolerance',        0.01);
options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfunccount @optimplotfval @optimplotstepsize @optimplotfirstorderopt });
%options = optimoptions(options,'ConstraintTolerance',  0.01);

for i = 1:SBOiter_max
    
    % copt
    x(i) = fmincon(@(x) norm(14+pdistr(x,U_mag,xp,yp,'delta')),[14],[],[],[],[],10,20,[],options);
       
    temp_fig = gcf;
    savefig(temp_fig,strcat('./Output/',num2str(i),'opt.fig'));
    close gcf
    
    fc{i} = pdistr(x(i),U_mag,xp,yp,'all');
    xc = fc{i}{3};
    ff{i} = xfoil_distr(x(i),xc,'all');
    
    figure(i);
    subplot(1,2,1)
    plot(xc,ff{i}{2},'go-')
    hold on
    plot(xc,fc{i}{2},'bo-')
    grid on;
    set(gca,'YDir','Reverse'); % reverse y axis
    
    legend('xfoil_v','HS')
    xlabel('corda');
    ylabel('-C_p');
    
    title(sprintf('alpha = %1.1f deg',x(i)));
    
    subplot(1,2,2)
    plot(xc,ff{i}{2}-fc{i}{2},'go-')
    grid on;
    xlabel('corda');
    ylabel('| DC_p |');
    title('f-c_{opt}')
    
    temp_fig = gcf;
    savefig(temp_fig,strcat('./Output/',num2str(i),'comp.fig'));
    close gcf
    
end

    
%% SUBROUTINE

function f = pdistr(alpha_degr,U_mag,xp,yp,wtd)

alpha = alpha_degr*pi/180; % incidenza profilo
U_inf = U_mag*[cos(alpha) sin(alpha)]; % velocità asintotica

[~,~,B,xc,~,~,~,T,x] = HS_staz(xp,yp,U_inf);
U_inf_t = U_inf*T;
u_pert = B*x + U_inf_t';

switch wtd
    case 'distro'
        f = (1-(u_pert./norm(U_inf)).^2)';
    case 'delta'
        cp = (1-(u_pert./norm(U_inf)).^2)';
        f  = min(cp - cp(1));  
    case 'all'
        cp = (1-(u_pert./norm(U_inf)).^2)';
        f{1}  = min(cp - cp(1));
        f{2}  = cp;
        f{3}  = xc;
    otherwise
        error('error')
end

end

function f = xfoil_distr(alpha_v,xc,wtd,case_number)

if nargin == 3
   case_number = 1;
end
[~,foil] = xfoil2matlab('NACA0012',alpha_v,9000000,0,1000);

n   = size(xc,2)/2;
ncp = size(foil.xcp,1)/2;

% reinterpolo
foil.cpI = fliplr([spline(foil.xcp(1:ncp),foil.cp(1:ncp),xc(1:n)),...;
    spline(foil.xcp(ncp+1:end),foil.cp(ncp+1:end),xc(n+1:end))]);

switch wtd
    case 'distro'
        f = foil.cpI;
    case 'delta'
        f = min(foil.cpI - foil.cpI(1));
    case 'all'
        f{1}  = min(foil.cpI - foil.cpI(1));
        f{2}  = foil.cpI;
    otherwise
        error('error')
end

% pulizia
d1 = system(sprintf...
    ('mv ./Routines/:00.bl ./Output/%d.bl',case_number));
d2 = system(sprintf...
    ('mv ./Routines/xfoil2matlab.inp ./Output/%d.inp',case_number));
d3 = system(sprintf...
    ('mv ./Routines/xfoil2matlab_pwrt.dat ./Output/%d.dat',case_number));
d4 = system(sprintf...
    ('mv ./Routines/xfoil.out ./Output/%d.out',case_number));

end

