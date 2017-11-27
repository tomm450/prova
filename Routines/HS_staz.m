%% HS
function [lfun,A,B,xc,yc,dl,N,T,x,theta_G] = HS_staz(xp,yp,U_inf)
%questo è giusto
lfun = localfunctions; %restituisce le funzioni locali come globali
if nargout == 1
    return
end

Npt = length(xp); % numero di punti
Npan = length(xp)-1; % numero di pannelli

%% Punti di collocazione
[xc,yc] = collocation_points(xp,yp);

%% Lunghezze pannelli
[dx,dy,dl] = panels_dimensions(xp,yp);

%% Versori normali e tangenti
[N,T] = NTversors(dx,dy,dl);

%% Angoli globali
[theta_G] = panels_inclination(xp,yp);

%% Velocità indotte
% Preallocazione matrici coefficienti
A = zeros(Npt,Npt);
C = A;
B = zeros(Npt-1,Npt);
D = B;
b = zeros(Npt,1);
for j = 1:Npan
    % Matrice di rotazione per portarci nel sistema pannello
    [R] = rotation_matrix(theta_G(j));
    for k = 1:Npan
        % Raggi fra estremi pannello j e punti di collocazione k
        % Angoli fra estremi pannello j e punti di collocazione k
        [r1,r2,dt] = radius_dtheta(xc(k),yc(k),xp(j),yp(j),xp(j+1),yp(j+1),R,j,k);
        % Velocità indotte SORGENTI
        [Us,Vs] = source_induced_speed(r1,r2,dt);
        
        % Velocità indotte VORTICI
        [Uv,Vv] = vortex_induced_speed(r1,r2,dt);
        
        % Vettori velocità GLOBALI
        S_vel = R'*[Us Vs]';
        V_vel = R'*[Uv Vv]';
        
        A(k,j) = S_vel'*N(:,k);
        C(k,j) = V_vel'*N(:,k);
        B(k,j) = S_vel'*T(:,k);
        D(k,j) = V_vel'*T(:,k);
        
    end
    % Vettore termini noti
    b(j,1) = -U_inf*N(:,j);
end

%% Somma velocità indotte dai vortici su ogni pannello
A(:,end) = sum(C,2);
B(:,end) = sum(D,2);

%% Condizione di Kutta
A(end,:) = B(1,:)+ B(end,:);
b(end,1) = -U_inf*(T(:,1)+T(:,end));

%% Risoluzione sistema lineare
if nargout == 8 %non è richiesta la risoluzione
    return
end

x = A\b;

end

%****************************************************************************%
% subfunction

function [xc,yc] = collocation_points(xp,yp)
xc = 0.5*(xp(2:end)+xp(1:end-1));
yc = 0.5*(yp(2:end)+yp(1:end-1));
end

function [dx,dy,dl] = panels_dimensions(xp,yp)
dx = xp(2:end)-xp(1:end-1);
dy = yp(2:end)-yp(1:end-1);
dl = sqrt(dx.^2+dy.^2);
end

function [N,T] = NTversors(dx,dy,dl)
N = [(-dy./dl);
    (dx./dl)];
N=N./sqrt(N(1,:).^2+N(2,:).^2);
T = [(dx./dl);
    (dy./dl)];
T=T./sqrt(T(1,:).^2+T(2,:).^2);
end

function [theta_G] = panels_inclination(xp,yp)
theta_G = atan2(yp(2:end)-yp(1:end-1),xp(2:end)-xp(1:end-1));
end

function [R] = rotation_matrix(theta)
R = [cos(theta)  sin(theta)
    -sin(theta) cos(theta)];
end

function [r1,r2,dt] = radius_dtheta(xc,yc,xp,yp,xp2,yp2,R,varargin)
if nargin == 7
    j=2;k=3;
else
    j = varargin{1};
    k = varargin{2};
end

r1_G = [xc-xp yc-yp];
r2_G = [xc-xp2 yc-yp2];
r1 = R*r1_G';
r2 = R*r2_G';

t1 = atan2(r1(2),r1(1));
t2 = atan2(r2(2),r2(1));
dt = t2-t1;
if j == k && t2-t1<1e-15
    dt=t1-t2;
end

end

function [Us,Vs] = source_induced_speed(r1,r2,dt,varargin)
if nargin == 3
    sigma = 1;
else
    sigma = varargin{1};
end

Us = -sigma./(2*pi)*log(norm(r2)/norm(r1));
Vs = sigma.*(dt)/(2*pi);
end

function [Uv,Vv] = vortex_induced_speed(r1,r2,dt,varargin)
if nargin == 3
    gamma = 1;
else
    gamma = varargin{1};
end

Uv = -gamma.*(dt)/(2*pi);
Vv = -gamma./(2*pi)*log(norm(r2)/norm(r1));
end



