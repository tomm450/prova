clc
clear all;
close all;

addpath('./Routines');
%% Definizione profilo

xfoil_cmd = '/home/tom/Downloads/Xfoil/bin/xfoil'; %portatile
%xfoil_cmd = 'xfoil';                               %fisso

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

SBOiter_max = 15;

% Start with the default options
OPT.nvars = 1;
options = optimoptions('fmincon');
% Modify options setting
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxFunctionEvaluations', 50);
options = optimoptions(options,'OptimalityTolerance', 0.005);
options = optimoptions(options,'FunctionTolerance',   0.005);
options = optimoptions(options,'StepTolerance',        0.01);
options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfunccount @optimplotfval @optimplotstepsize @optimplotfirstorderopt });
%options = optimoptions(options,'ConstraintTolerance',  0.01);

fvals_c = {};
fvals_f = {};
Skf     = {};

x = nan(1,SBOiter_max);
fcopt = cell(SBOiter_max,1);
fc = fcopt;
ff = fc;

flap_flag    = 1;
flap_perc    = 0.8; % hinge position (normalized respect to chord)
flap_deflect = 10;  % deflection degree

str_flap = sprintf('gdes flap %1.2f 0 %d exec', flap_perc, flap_deflect);

for iter_f = 1:SBOiter_max
    fprintf('\n\n\n#####################################################\n');
    % copt
    x(iter_f) = fmincon(@(x) norm(14+pdistr(x,U_mag,xp,yp,'delta',fvals_c,fvals_f,Skf)),...
        14,[],[],[],[],10,20,[],options);
    
    temp_fig = gcf;
    savefig(temp_fig,strcat('./Output/',num2str(iter_f),'opt.fig'));
    close gcf
    
    
    fc{iter_f}      = pdistr(x(iter_f),U_mag,xp,yp,'all');
    fcopt{iter_f}  = pdistr(x(iter_f),U_mag,xp,yp,'all',fvals_c,fvals_f,Skf);
    
    xc = fc{iter_f}{3};
   
    ff{iter_f} = xfoil_distr(x(iter_f),xc,'all',xfoil_cmd,iter_f,opt_str,str_flap);
    
    Efc(iter_f,:)    =  ff{iter_f}{2}-fc{iter_f}{2};
    Efcopt(iter_f,:) =  ff{iter_f}{2}-fcopt{iter_f}{2};
    
    fprintf('Errore medio F-C    = %f \n',mean(ff{iter_f}{2}-fc{iter_f}{2}));
    fprintf('Errore medio F-Copt = %f \n',mean(ff{iter_f}{2}-fcopt{iter_f}{2}));
    
    figure(iter_f);
    
    subplot(1,2,1)
    plot(xc,ff{iter_f}{2},'bo-')
    hold on
    plot(xc,fcopt{iter_f}{2},'ro-')
    plot(xc,fc{iter_f}{2},'go-')
    
    grid on;
    set(gca,'YDir','Reverse'); % reverse y axis
    
    legend('xfoil','HS_{cor}','HS')
    xlabel('corda');
    ylabel('-C_p');
    
    title(sprintf('alpha = %1.1f deg',x(iter_f)));
    
    subplot(1,2,2)
    plot(xc,ff{iter_f}{2}-fc{iter_f}{2},'bo-')
    hold on
    plot(xc,ff{iter_f}{2}-fcopt{iter_f}{2},'ro-')
    legend('xfoil - HS','xfoil - HS_{cor}')
    
    
    grid on;
    xlabel('corda');
    ylabel('| DC_p |');
    title('f-c_{opt}')
    
    temp_fig = gcf;
    savefig(temp_fig,strcat('./Output/',num2str(iter_f),'comp.fig'));
    close gcf
    
    fvals_copt{iter_f} = fcopt{iter_f}{2};
    fvals_c{iter_f}    = fc{iter_f}{2};
    fvals_f{iter_f}    = ff{iter_f}{2};
    
        %% correzione
    if iter_f == 1  % ho calcolato un fvals_c e un fvals_f, Skk e Skf saranno eye
        Skf = {eye(max(size(fvals_f{end})))};     % ho solo 1 valore di fvals, la prima CORREZIONE
        % Skh = {1};
        
        DF = [];
        DC = [];
        
    else
        clear DF DC
               
        %for w = 1: max(1,min([OPT.nvars-1,(size(fvals_f,2) -1)]))
        for w = 1: size(fvals_f,2) -1
             
            DF(:,w) = [fvals_f{end} - fvals_f{end-w}];
            DC(:,w) = [fvals_c{end} - fvals_c{end-w}];
            %             DISf(:,w) = [ISsol_f{end}(1) - ISsol_f{end-w}(1)];
            %             DISc(:,w) = [ISsol_c{end}(1) - ISsol_c{end-w}(1)];
            %
            %             DHf(:,w)  = [Hf{end} - Hf{end-w}];
            %             DHc(:,w)  = [Hc{end} - Hc{end-w}];
            
            
        end
        
        fprintf('DF(:,w) = [fvals_f{end} - fvals_f{end-w}];\nDC(:,w) = [fvals_c{end} - fvals_c{end-w}];\n')
        
        %DC
        
        %DF
        
        [Uf,Sf,Vf] = svd(DC);
        %         [Us,Ss,Vs] = svd(DISc);
        %         [Uh,Sh,Vh] = svd(DHc);
        
        
        % pseudoinversa
        Sf_cros = Sf';
        
        for j = 1: min([size(Sf_cros,1),size(Sf_cros,2)])
            if abs(Sf_cros(j,j)) >= 1e-10
                Sf_cros(j,j) = 1/Sf_cros(j,j);
            end
        end
        
        %         Ss_cros = Ss';
        %         for j = 1: min([size(Ss_cros,1),size(Ss_cros,2)])
        %             if abs(Ss_cros(j,j)) >= 1e-10
        %                 Ss_cros(j,j) = 1/Ss_cros(j,j);
        %             end
        %         end
        %
        %         Sh_cros = Sh';
        %         for j = 1: min([size(Sh_cros,1),size(Sh_cros,2)])
        %             if abs(Sh_cros(j,j)) >= 1e-10
        %                 Sh_cros(j,j) = 1/Sh_cros(j,j);
        %             end
        %         end
              
        DCc   = Vf*Sf_cros*Uf';
        %         DIScc = Vs*Ss_cros*Us';
        %         DHcc = Vh*Sh_cros*Uh';
        
        %T   = DC*DFc;
        
        Skf{end+1} = DF*DCc;
        %         Skk{end+1} = DISf*DIScc;
        %         Skh{end+1} = DHf*DHcc;
        
        
        %         Skf{end}
        %         Skk{end}
        %         Skh{end}
        
    end
    
    
end

% colorbar
n_col   = ceil(2+size(fvals_f,2)/2);

n_shade = linspace(0,1,n_col); 
n_shade = n_shade(2:end-1); 
n_shade = n_shade';

n_zeros = zeros(size(n_shade));
CMAT = [[1              ,0              ,      0]; ...
        [flipud(n_shade),n_shade        ,n_zeros]; ...
        [0              ,1              ,      0]; ...
        [n_zeros        ,flipud(n_shade),n_shade]];


for k = 1:size(fvals_f,2)
    f100 = figure(100);
    
    set(f100,'Position',[50 50 1200 800]);
    subplot(2,2,1)
    plot(Efcopt(k,:),'Color',CMAT(k,:)); 
    hold on;
    grid on;
    title('Evolution err(xfoil - HS_{cor})');   
    
    subplot(2,2,3)
    plot(Efc(k,:),'Color',CMAT(k,:)); 
    hold on;
    grid on;
    title('Evolution err(xfoil - HS)');
    
    subplot(2,2,2)
    semilogy(abs(Efcopt(k,:)),'Color',CMAT(k,:)); 
    hold on;
    grid on;
    title('Evolution log(xfoil - HS_{cor})');
    
    subplot(2,2,4)
    semilogy(abs(Efc(k,:)),'Color',CMAT(k,:)); 
    hold on;
    grid on;
    title('Evolution log(xfoil - HS)');
    
    leg_cell{k} = num2str(k);
    
    
end

figure(100);
% subplot(2,2,1)
% legend(leg_cell)
subplot(2,2,3)
legend(leg_cell)
% subplot(2,2,3)
% legend(leg_cell)
% subplot(2,2,4)
% legend(leg_cell)

savefig(f100,'./Output/evo.fig');

f101 = figure(101);
plot(x,'o--');
grid on
title('Alpha opt');
xlabel('Iteration');
ylabel('Alpha [deg]');
savefig(f101,'./Output/alpha_conv.fig');

%% SUBROUTINE

function f = pdistr(alpha_degr,U_mag,xp,yp,wtd,fvals_c,fvals_f,Skf)

if nargin == 5
    fvals_c  = {};
    fvals_f  = {};
    Skf      = {};
end

alpha = alpha_degr*pi/180; % incidenza profilo
U_inf = U_mag*[cos(alpha) sin(alpha)]; % velocità asintotica

[~,~,B,xc,~,~,~,T,x] = HS_staz(xp,yp,U_inf);
U_inf_t = U_inf*T;
u_pert = B*x + U_inf_t';

switch wtd
    case 'distro'
                       
        f = (1-(u_pert./norm(U_inf)).^2)';
        
        if size(fvals_f,2) > 0
           f = fvals_f{end}' + Skf{end}*( f' - fvals_c{end}' );    
           f = f';
        end
        
        
    case 'delta'
    
        cp = (1-(u_pert./norm(U_inf)).^2)';
        
        if size(fvals_f,2) > 0
           cp = fvals_f{end}' + Skf{end}*( cp' - fvals_c{end}' );    
           cp = cp';
        end
        
        f  = min(cp - cp(1));  
    
    case 'all'
        
        cp = (1-(u_pert./norm(U_inf)).^2)';
        
        if size(fvals_f,2) > 0
           cp = fvals_f{end}' + Skf{end}*( cp' - fvals_c{end}' );    
           cp = cp';
        end
        
        f{1}  = min(cp - cp(1));
        f{2}  = cp;
        f{3}  = xc;
    
    otherwise
        error('error')
end

end

function f = xfoil_distr(alpha_v,xc,wtd,xfoil_cmd,case_number,opt_str)

if nargin == 4
   case_number = 1;
   opt_str = [];
end

if size(opt_str,2) == 0
    [~,foil] = xfoil2matlab('NACA0012',alpha_v,9000000,0,1000,xfoil_cmd);
else
    [~,foil] = xfoil2matlab('NACA0012',alpha_v,9000000,0,1000,xfoil_cmd,opt_str);
end
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

