clc
clear all; 
close all;
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

alpha_v = [12:14];

U_mag   = 135;

iter_number = 5000;
 
cp  = nan(size(alpha_v,2),Npan);
DCpHS    = nan(size(alpha_v,2),1); 
DCpFoil  = nan(size(alpha_v,2),1); 
DCpFoilv  = nan(size(alpha_v,2),1); 
DCpFoilo  = nan(size(alpha_v,2),1); 

for av = 1:size(alpha_v,2)

fprintf('alpha = %1.2f deg \n\n',alpha_v(av));

alpha = alpha_v(av)*pi/180; % incidenza profilo
U_inf = U_mag*[cos(alpha) sin(alpha)]; % velocità asintotica

[~,A,B,xc,yc,dl,N,T,x] = HS_staz(xp,yp,U_inf);
U_inf_t = U_inf*T;
u_pert = B*x + U_inf_t';

%% Coefficiente di pressione
cp(av,:) = (1-(u_pert./norm(U_inf)).^2)';
DCpHS(av)  = min(cp(av,:) - cp(av,1));

end

% [pol,foil] = xfoil2matlab(coord,alpha,Re,Mach,varargin)
[polo,foilo] = xfoil2matlab('NACA0012',alpha_v,9000000,0,1000);
[pol,foil]   = xfoil2matlab([xc;yc],alpha_v,0,0,5000);
[polv,foilv] = xfoil2matlab([xc;yc],alpha_v,9000000,0,1000);



for k = 1:size(alpha_v,2)
    
    n   = size(xc,2)/2;
    ncp = size(foilv.xcp,1)/2;
    % reinterpolo
    foilv.cpI(:,k) = fliplr([spline(foilv.xcp(1:ncp),foilv.cp(1:ncp,k),xc(1:n)),...;
                     spline(foilv.xcp(ncp+1:end),foilv.cp(ncp+1:end,k),xc(n+1:end))]);
    
    foilo.cpI(:,k) = fliplr([spline(foilo.xcp(1:ncp),foilo.cp(1:ncp,k),xc(1:n)),...;
                     spline(foilo.xcp(ncp+1:end),foilo.cp(ncp+1:end,k),xc(n+1:end))]);
    
    foil.cpI(:,k) = fliplr([spline(foil.xcp(1:ncp),foil.cp(1:ncp,k),xc(1:n)),...;
                     spline(foil.xcp(ncp+1:end),foil.cp(ncp+1:end,k),xc(n+1:end))]);
                 
    
    Emat{k} = [cp(k,:)'-foilv.cpI(:,k) ,cp(k,:)'-foilo.cpI(:,k),cp(k,:)'-foil.cpI(:,k)];           
                 
    figure(k+1);
    subplot(1,2,1)
    plot(xc,foilv.cpI(:,k),'go-')
    hold on
    plot(xc,foilo.cpI(:,k),'co-')
    plot(xc,foil.cpI(:,k),'bo-')
    plot(xc,cp(k,:),'ro-')
    grid on;
    set(gca,'YDir','Reverse'); % reverse y axis
         
    legend('xfoil_v','native coordinates','xfoil_i','HS')
    xlabel('corda');
    ylabel('-C_p');
    
    title(sprintf('alpha = %1.1f deg',alpha_v(k)));
    
    subplot(1,2,2)
    semilogy(xc,abs(Emat{k}(:,1)),'go-')
    hold on
    semilogy(xc,abs(Emat{k}(:,2)),'co-')
    semilogy(xc,abs(Emat{k}(:,3)),'bo-')
    
    grid on;
    alfa=num2str(alpha_v(av)*180/pi);
    legend('xfoil_v-HS','native coordinates-HS','xfoil_i-HS')
    xlabel('corda');
    ylabel('| DC_p |');
    
    tit1 = 'COEFFICIENTE DI PRESSIONE';
    tit2 = 'Differenza modelli';
    title({tit1;tit2})
    
    DCpFoil(k)  = min(foil.cp(:,k) - foil.cp(1,k));
    DCpFoilv(k) = min(foilv.cp(:,k) - foilv.cp(1,k));
    DCpFoilo(k) = min(foilo.cp(:,k) - foilo.cp(1,k));
    
end

