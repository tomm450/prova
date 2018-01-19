clc
clear all;
close all;

addpath('./Routines');
%% Definizione profilo
naca = 'NACA0012';
global REYNOLDS
REYNOLDS = 9000000

n = 160; % n pannelli 
[xp,yp] = NACA_generator(naca,n,'cos',0);

% naca = load('64212.mat');
% dummy = size(naca.air,1); dummy = [1,ceil(linspace(2,dummy-1,n*2-3)),dummy];
% xp = naca.air(dummy,1)';
% yp = naca.air(dummy,2)';

%native = load('./nativeNACA0012.txt');

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

xfoil_cmd = '/home/tom/Downloads/Xfoil/bin/xfoil'; %portatile
%xfoil_cmd = 'xfoil';                               %fisso


%% cose normalmente non note
alpha_test = [0:30];
% 
% % polare Xfoil inv
% [polI,foilI] = xfoil2matlab(naca,alpha_test,0,0,5000,xfoil_cmd);
% polI.Eff = polI.CL./polI.CD;
% % polare Xfoil visc
% [pol,foil] = xfoil2matlab(naca,alpha_test,REYNOLDS,0,5000,xfoil_cmd);
% pol.Eff = pol.CL./pol.CD;
% 
% 
% for i = 1:size(alpha_test,2)
%    % calcolo deltaP xfoil 
%    pol.DeltaP(i)  = min(foil.cp(:,i))-foil.cp(end,i); 
%    polI.DeltaP(i) = min(foilI.cp(:,i))-foilI.cp(end,i); 
%    % calcolo polare HS
%    fprintf('HS alpha = %d deg; \n',alpha_test(i));
%    f = pdistr(alpha_test(i),U_mag,xp,yp,'all');
%    hs.CL(i) = f{4}(1);
%    hs.CD(i) = f{4}(2);
%    hs.Eff(i) = f{4}(1)/f{4}(2);
%    hs.DeltaP(i) = f{1};
%    
% end

load('HSdata0012.mat')
load('FOILdata0012.mat')

[CLmax,iCL] = max(pol.CL);

f1 = figure(1234);

subplot(2,2,1)
plot(pol.alpha,pol.CL,'b',polI.alpha,polI.CL,'c',pol.alpha,hs.CL,'r')
title('Cl'); grid on

subplot(2,2,2)
yyaxis left;
plot(pol.alpha,pol.CD,'b')
yyaxis right;
plot(pol.alpha,hs.CD,'r'); %,polI.alpha,polI.CD,'c');
title('Cd'); grid on

subplot(2,2,3)
yyaxis left;
plot(pol.alpha,pol.Eff,'b');
yyaxis right;
plot(pol.alpha,hs.Eff,'r'); %,polI.alpha,polI.Eff,'c'
title('Efficieza'); grid on
%axis([-1 30 0 150])

subplot(2,2,4)
plot(pol.alpha,abs(pol.DeltaP),'b',polI.alpha,abs(polI.DeltaP),'c',pol.alpha,abs(hs.DeltaP),'r');
hold on
plot(pol.alpha(iCL),abs(pol.DeltaP(iCL)),'bo',...
     polI.alpha(iCL),abs(polI.DeltaP(iCL)),'co',...
     pol.alpha(iCL),abs(hs.DeltaP(iCL)),'ro')
hold on
legend('Xfoil','Xfoil inv','HS',...
    sprintf('Xfoil|_{maxCl}     = %2.3f',abs(pol.DeltaP(iCL))),...
    sprintf('Xfoil inv|_{maxCl} = %2.3f',abs(polI.DeltaP(iCL))),...
    sprintf('HS|_{maxCl}        = %2.3f',abs(hs.DeltaP(iCL))),'Location','best')
plot(pol.alpha,14*ones(size(pol.alpha)),'m--')
set(f1,'Position',[10 10 1200 700])
title('abs(deltaP)')
grid on

SBOiter_max = 10;

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
x0 = [5:5:15];
x     = nan(1,SBOiter_max);
fcopt = cell(SBOiter_max,1);
fc    = cell(SBOiter_max,1);
ff    = cell(SBOiter_max,1);

Merr = 1;
iter_f = 1;
%for iter_f = 1:SBOiter_max

while iter_f-1 < SBOiter_max && abs(Merr) > 1e-5
    
    fprintf('\n\n\n#####################################################\n');
    fprintf('%d\n',iter_f);
    % copt
    %x(iter_f) = fmincon(@(x) norm(14+pdistr(x,U_mag,xp,yp,'delta',fvals_c,fvals_f,Skf)),...
    %    14,[],[],[],[],10,20,[],options);
    
    fitness = @(x) -pdistr(x,U_mag,xp,yp,'cl',fvals_c,fvals_f,Skf);%...
%                  ./pdistr(x,U_mag,xp,yp,'cd',fvals_c,fvals_f,Skf);
                  
    nonlinearcostr = @(x) deal(-14-pdistr(x,U_mag,xp,yp,'delta',fvals_c,fvals_f,Skf),0); 
    
    % chiamo fmincon
    for i = 1:size(x0,2)
    fprintf('OPT PARTENDO DA %1.2f (caso %d/%d) \n',x0(i),i,size(x0,2));
    [x_ms(iter_f,i),f_ms(iter_f,i)] = fmincon(fitness,...
        x0(i),[],[],[],[],0,20,...
        nonlinearcostr,...
        options);
    
    temp_fig = gcf;
    savefig(temp_fig,strcat('./Output/',num2str(iter_f),'opt',num2str(i),'.fig'));
    close gcf   
    end
    
    
    [f_ms_min,i_min] = min(f_ms(iter_f,:));
    
    x(iter_f) = x_ms(iter_f,i_min);
    % high res
    alpha_double_check = 2*linspace(floor(x(iter_f)/2),ceil(x(iter_f)/2),6);
    
    for i = 1:size(alpha_test,2)

        % calcolo polare HS
        fprintf('HS cor alpha = %d deg; (caso %d/%d)\n',alpha_test(i),i,size(alpha_test,2));
        f = pdistr(alpha_test(i),U_mag,xp,yp,'all',fvals_c,fvals_f,Skf);
        surPol.CL{iter_f}(i) = f{4}(1);
        surPol.CD{iter_f}(i) = f{4}(2);
        surPol.Eff{iter_f}(i) = f{4}(1)/f{4}(2);
        surPol.DeltaP{iter_f}(i) = f{1};
        
    end
    
    for i = 1:size(alpha_double_check,2)

        % calcolo polare HS
        fprintf('HS double check alpha = %d deg; (caso %d/%d)\n',alpha_double_check(i),i,size(alpha_double_check,2));
        f = pdistr(alpha_test(i),U_mag,xp,yp,'all',fvals_c,fvals_f,Skf);
        double_check.CL{iter_f}(i) = f{4}(1);
        double_check.CD{iter_f}(i) = f{4}(2);
        double_check.Eff{iter_f}(i) = f{4}(1)/f{4}(2);
        double_check.DeltaP{iter_f}(i) = f{1};
        
    end
    double_check.alpha{iter_f}  = alpha_double_check;
    
    
  
    
    
    fc{iter_f}     = pdistr(x(iter_f),U_mag,xp,yp,'all');
    fcopt{iter_f}  = pdistr(x(iter_f),U_mag,xp,yp,'all',fvals_c,fvals_f,Skf);   
        
    xc = fc{iter_f}{3};
    
    ff{iter_f} = xfoil_distr(x(iter_f),xc,'all',xfoil_cmd,iter_f);
    
    
    
    fprintf('Confronto onboard (c,copt,f)...\n');
    fprintf('deltaP = %+02.3f, %+02.3f, %+02.3f \n',fc{iter_f}{1},fcopt{iter_f}{1},ff{iter_f}{1});
    fprintf('Cl     = %+02.3f, %+02.3f, %+02.3f \n',fc{iter_f}{end}(1),fcopt{iter_f}{end}(1),ff{iter_f}{end}(1));
    fprintf('Cd     = %+02.3f, %+02.3f, %+02.3f \n',fc{iter_f}{end}(2),fcopt{iter_f}{end}(2),ff{iter_f}{end}(2));
    
    [nlci,nlce] = nonlinearcostr( x(iter_f));
    fprintf('Vincolo non lineare -> %f \n\n',nlci);    
    Efc(iter_f,:)    =  ff{iter_f}{2}-fc{iter_f}{2};
    Efcopt(iter_f,:) =  ff{iter_f}{2}-fcopt{iter_f}{2};
    
    fprintf('Errore medio F-C    = %f \n',mean(ff{iter_f}{2}-fc{iter_f}{2}));
    fprintf('Errore medio F-Copt = %f \n',mean(ff{iter_f}{2}-fcopt{iter_f}{2}));
    
    Merr = mean(ff{iter_f}{2}-fcopt{iter_f}{2});
    
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
    iter_f = iter_f +1;
    
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

save('maxEff.mat');
dummy = 22



%% SUBROUTINE

function f = pdistr(alpha_degr,U_mag,xp,yp,wtd,fvals_c,fvals_f,Skf)

if nargin == 5
    fvals_c  = {};
    fvals_f  = {};
    Skf      = {};
end

alpha = alpha_degr*pi/180; % incidenza profilo
U_inf = U_mag*[cos(alpha) sin(alpha)]; % velocità asintotica

[~,~,B,xc,~,dl,~,T,x,theta_G] = HS_staz(xp,yp,U_inf);
U_inf_t = U_inf*T;
u_pert  = B*x + U_inf_t';

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
        
    case 'cl'
        
        cp = (1-(u_pert./norm(U_inf)).^2)';
        
        if size(fvals_f,2) > 0
            cp = fvals_f{end}' + Skf{end}*( cp' - fvals_c{end}' );
            cp = cp';
        end
        
        n = size(xc,2)/2;
        
        cp_air  = cp(1:2*n);
        %cp_slat = cp(2*n+1:end);
        
        CD_ref_air = sum(  dl(1,:).*cp_air.*sin(theta_G(1,:)));% + dl(2,:).*cp_slat'.*sin(theta_G(2,:)));
        CL_ref_air = sum( -dl(1,:).*cp_air.*cos(theta_G(1,:)));% - dl(2,:).*cp_slat'.*cos(theta_G(2,:)));
        
        CL = CL_ref_air*cosd(alpha_degr) - CD_ref_air*sind(alpha_degr);
        %CD = CL_ref_air*sind(alpha_degr) + CD_ref_air*cosd(alpha_degr);
        
        if size(CL) ~= [1,1]
            error('dimensioni CL errate!');
        end
        %         if size(CD) ~= [1,1]
        %             error('dimensioni CD errate!');
        %         end
        f = [CL];
        
    case 'cd'
        
        cp = (1-(u_pert./norm(U_inf)).^2)';
        
        if size(fvals_f,2) > 0
            cp = fvals_f{end}' + Skf{end}*( cp' - fvals_c{end}' );
            cp = cp';
        end
        
        n = size(xc,2)/2;
        
        cp_air  = cp(1:2*n);
        %cp_slat = cp(2*n+1:end);
        
        CD_ref_air = sum(  dl(1,:).*cp_air.*sin(theta_G(1,:)));% + dl(2,:).*cp_slat'.*sin(theta_G(2,:)));
        CL_ref_air = sum( -dl(1,:).*cp_air.*cos(theta_G(1,:)));% - dl(2,:).*cp_slat'.*cos(theta_G(2,:)));
        
        %CL = CL_ref_air*cosd(alpha_degr) - CD_ref_air*sind(alpha_degr);
        CD = CL_ref_air*sind(alpha_degr) + CD_ref_air*cosd(alpha_degr);
        
        %         if size(CL) ~= [1,1]
        %             error('dimensioni CL errate!');
        %         end
        if size(CD) ~= [1,1]
            error('dimensioni CD errate!');
        end
        f = [CD];
        
    case 'coeff'
        
        cp = (1-(u_pert./norm(U_inf)).^2)';
        
        if size(fvals_f,2) > 0
            cp = fvals_f{end}' + Skf{end}*( cp' - fvals_c{end}' );
            cp = cp';
        end
        
        n = size(xc,2)/2;
        
        cp_air  = cp(1:2*n);
        %cp_slat = cp(2*n+1:end);
        
        CD_ref_air = sum(  dl(1,:).*cp_air.*sin(theta_G(1,:)));% + dl(2,:).*cp_slat'.*sin(theta_G(2,:)));
        CL_ref_air = sum( -dl(1,:).*cp_air.*cos(theta_G(1,:)));% - dl(2,:).*cp_slat'.*cos(theta_G(2,:)));
        
        CL = CL_ref_air*cosd(alpha_degr) - CD_ref_air*sind(alpha_degr);
        CD = CL_ref_air*sind(alpha_degr) + CD_ref_air*cosd(alpha_degr);
        
        if size(CL) ~= [1,1]
            error('dimensioni CL errate!');
        end
        if size(CD) ~= [1,1]
            error('dimensioni CD errate!');
        end
        f = [CL,CD];
        
    case 'all'
        
        cp = (1-(u_pert./norm(U_inf)).^2)';
        
        if size(fvals_f,2) > 0
            cp = fvals_f{end}' + Skf{end}*( cp' - fvals_c{end}' );
            cp = cp';
        end
        
        n = size(xc,2)/2;
        
        cp_air  = cp(1:2*n);
        %cp_slat = cp(2*n+1:end);
        
        CD_ref_air = sum(  dl(1,:).*cp_air.*sin(theta_G(1,:)));% + dl(2,:).*cp_slat'.*sin(theta_G(2,:)));
        CL_ref_air = sum( -dl(1,:).*cp_air.*cos(theta_G(1,:)));% - dl(2,:).*cp_slat'.*cos(theta_G(2,:)));
        
        CL = CL_ref_air*cosd(alpha_degr) - CD_ref_air*sind(alpha_degr);
        CD = CL_ref_air*sind(alpha_degr) + CD_ref_air*cosd(alpha_degr);
        
        if size(CL) ~= [1,1]
            error('dimensioni CL errate!');
        end
        if size(CD) ~= [1,1]
            error('dimensioni CD errate!');
        end
        
        f{1}  = min(cp - cp(1));
        f{2}  = cp;
        f{3}  = xc;
        f{4} = [CL,CD];
        
    otherwise
        error('error')
end

end

function f = xfoil_distr(alpha_v,xc,wtd,xfoil_cmd,case_number)

global REYNOLDS

if nargin == 4
   case_number = 1;
end

%[pol,foil] = xfoil2matlab('NACA0012',alpha_v,9000000,0,1000,xfoil_cmd);
[pol,foil] = xfoil2matlab('NACA0012',alpha_v,REYNOLDS,0,1000,xfoil_cmd);


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
    case 'cl'
        f = pol.CL;
    case 'cd'
        f = pol.CD;
    case 'coeff'
        f = [pol.CL,pol.CD];
    case 'all'
        f{1}  = min(foil.cpI - foil.cpI(1));
        f{2}  = foil.cpI;
        f{3}  = [pol.CL,pol.CD];
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

