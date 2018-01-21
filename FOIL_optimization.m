clc
clear all;
close all;

addpath('./Routines');
%% Definizione profilo
naca = 'NACA0012';
global REYNOLDS
REYNOLDS = 9000000


iter_number = 5000;

xfoil_cmd = '/home/tom/Downloads/Xfoil/bin/xfoil'; %portatile
%xfoil_cmd = 'xfoil';                               %fisso


% Start with the default options
OPT.nvars = 1;
options = optimoptions('fmincon');
% Modify options setting
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxFunctionEvaluations', 50);
options = optimoptions(options,'OptimalityTolerance', 0.005);
options = optimoptions(options,'FunctionTolerance',   0.005);
options = optimoptions(options,'StepTolerance',        0.001);
options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfunccount @optimplotfval @optimplotstepsize @optimplotfirstorderopt });
%options = optimoptions(options,'ConstraintTolerance',  0.01);

x0 = [5:5:15];

    
    fprintf('\n\n\n#####################################################\n');
    
    fitness = @(x) -xfoil_distr(x,'cl',xfoil_cmd,999)...
                  ./xfoil_distr(x,'cd',xfoil_cmd,999);
    a = [10:0.05:12]; 
    E = [];
    A = [];
    for k = 1:max(size(a))
        e = fitness(a(k));
        fprintf('a = %2.2f e = %4.4f \n',a(k),e);
        if size(e) == [0 0]
        else
            E = [E;e];
            A = [A;a(k)];
        end
    end
    
%     nonlinearcostr = @(x) deal(-14-xfoil_distr(x,'delta',xfoil_cmd,999),0); 
%     
%     % chiamo fmincon
%     for i = 1:size(x0,2)
%     fprintf('OPT PARTENDO DA %1.2f (caso %d/%d) \n',x0(i),i,size(x0,2));
%     [x_ms(1,i),f_ms(1,i)] = fmincon(fitness,...
%         x0(i),[],[],[],[],0,20,...
%         nonlinearcostr,...
%         options);
%     
%     temp_fig = gcf;
%     savefig(temp_fig,strcat('./Output/',num2str(999),'opt',num2str(i),'.fig'));
%     close gcf   
%     end
%     
%     
%     [f_ms_min,i_min] = min(f_ms(1,:));
%     
%     x(1) = x_ms(iter_f,i_min);
%     

%% SUBROUTINE

function f = xfoil_distr(alpha_v,wtd,xfoil_cmd,case_number)

global REYNOLDS

if nargin == 4
   case_number = 1;
end

%[pol,foil] = xfoil2matlab('NACA0012',alpha_v,9000000,0,1000,xfoil_cmd);
[pol,foil] = xfoil2matlab('NACA0012',alpha_v,REYNOLDS,0,1000,xfoil_cmd);


% n   = size(xc,2)/2;
% ncp = size(foil.xcp,1)/2;
% 
% % reinterpolo
% foil.cpI = fliplr([spline(foil.xcp(1:ncp),foil.cp(1:ncp),xc(1:n)),...;
%     spline(foil.xcp(ncp+1:end),foil.cp(ncp+1:end),xc(n+1:end))]);

switch wtd
    case 'distro'
        f = foil.cp;
    case 'delta'
        f = min(foil.cp - foil.cp(1));
    case 'cl'
        f = pol.CL;
    case 'cd'
        f = pol.CD;
    case 'coeff'
        f = [pol.CL,pol.CD];
    case 'all'
        f{1}  = min(foil.cp - foil.cp(1));
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

