function [xx,yy] = NACA_generator(NACA,N,varargin )
%NACA_Airfoil: Obtains the characteristic for a NACA 4-digit and 5-digit 
%               series and plot the shape.
%               
%           Ex. [xx,yy] = NACA_generator( NACA,N,option)
%               
%               Where the first inpute is the NACA followed by a 4 or 5 
%               number digit as a string, NACAXXXX. The second input the number of 
%               points to plot the airfoil. The options are optional 
% 
%               first option  | interp          | 'linear' (default) or
%               'cos' 
%               second option | plot | 1 (default), 0
% 
% Ex:
%   NACA_Airfoil( 'NACA0012',30,'cos',1 )
% is a NACA 0012 evaluated in 30 points along x with cos method. The
% airfoil is plotted
% 
% ------------------------------------------------------------------------
% Author: Andrea Fani
% Date created: 08/12/2016

% 

%% Error
L = length(NACA);

%% Extract data
if L == 8
% 4-digit NACA airfoil
m = str2double(NACA(5))/100; % maximun chamber
p = str2double(NACA(6))/10; % maximun chamber along the cord
t = str2double(NACA(7:8))/100; % max thickness
elseif L == 9
% 5-digit NACA airfoil
cl = str2double(NACA(5))*3/20; % design coefficeint of lift
p = str2double(NACA(6))/20; % position of maximun chamber
q = str2double(NACA(7)); % check if it is a reflex camber
t = str2double(NACA(8:9))/100; % max thickness
end

% Checks for trailing edge close or open and if it wants to be ploted
if ~isempty(varargin)
    xint = varargin{1};                 % varagin non vuoto
    if length(varargin) == 2            
        
        PLT = varargin{2};              % due input, flap_info non esite
        flap_info = [];
    
    elseif length(varargin) == 3    
        PLT = varargin{2};
        flap_info = varargin{3};        % tre input, ho tutto
    
    elseif length(varargin) > 3         % errore negli input
        error('Too many input arguments')
    
    else
        PLT = 1;
        flap_info = [];    
    end
else
    xint='linear';
    PLT = 1;
    flap_info = [];
end

switch lower(xint)
    case 'linear'
        x=linspace(0,1,N);
    case 'cos'
        x=0.5*(1-cos(pi*linspace(0,1,N)));
    otherwise 
        error('xint must be linear or cos')
end

if size(flap_info,2) > 0
   
    if sum(x==flap_info(1)) > 0 
       % ho già il punto nel punto giusto
       [~,imin] = min(abs(x-flap_info(1)));
    else
        [~,imin] = min(abs(x-flap_info(1)));
        if imin == 1 || imin == size(x,2)
            error('aggiungere punti o modificare flap_info')
        end
        x = [x(1:imin-1),flap_info(1),x(imin+1:end)];
    end
end
%% Shape of mean camber
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1036; % if close


yt = (t/0.2).*(a0.*sqrt(x) + a1.*(x) + a2.*(x).^2 +...
    a3.*(x).^3 + a4.*(x).^4);

%%  camber line
if L == 8
    [yc,zeta] = camber('4digit',p,x,N,m);
elseif L == 9
    [yc,zeta] = camber('5digit',p,x,N,q,cl);
end

x_U = (x - yt.*sin(zeta));
x_L = (x + yt.*sin(zeta));

y_U = (yc + yt.*cos(zeta));
y_L = (yc - yt.*cos(zeta));

if size(flap_info,2)>0
    
    % porzione "buona"
    x_Ut = x_U(1:imin);
    x_Lt = x_L(1:imin);
    y_Ut = y_U(1:imin);
    y_Lt = y_L(1:imin);
        
    R = [ cosd(-flap_info(2)) -sind(-flap_info(2));...
          sind(-flap_info(2))  cosd(-flap_info(2))];
    
    RU = R*[x_U(imin+1:end)-flap_info(1);...
        y_U(imin+1:end)];
    
    RL = R*[x_L(imin+1:end)-flap_info(1);...
        y_L(imin+1:end)];
         
    
    jumpO = x(imin) - x(imin-1);
    jumpL = flap_info(1)+RL(1,1) - x_Lt(end);
    trim = 0;
    
    while jumpL <= 0.5*jumpO
       trim = trim+1;
       x_Lt = x_Lt(1:end-trim);
       y_Lt = y_Lt(1:end-trim);
       
       x_Ut = x_Ut(1:end-trim);
       y_Ut = y_Ut(1:end-trim);
       
       jumpL = flap_info(1)+RL(1,1) - x_Lt(end);
     
    end
    
    x_U = [x_Ut,flap_info(1)+RU(1,:)];
    y_U = [y_Ut,RU(2,:)];
    
    x_L = [x_Lt,flap_info(1)+RL(1,:)];
    y_L = [y_Lt,RL(2,:)];
end
  
xx=[x_L(end:-1:1),x_U(2:end)];
yy=[y_L(end:-1:1),y_U(2:end)];

%% Plot
if PLT == 1
figure;
plot(x_U,y_U,'bx-','LineWidth',2); hold on
plot(x_L,y_L,'bx-','LineWidth',2);
title([NACA(1:4) ' ' NACA(5:end)])

% domain
% xmin = -.1;
% xmax = c*.1+c;
% ymin = -2*max(Upper(:,2));
% ymax = 2*max(Upper(:,2));

% axis([xmin xmax ymin ymax])
axis equal
end
end

%% Compute the camber for different airfoils
function [yc,zeta] = camber(naca,p,x,N,varargin)
I = find(x > p,1);
yc = zeros(1,N);
dyc_dx = zeros(1,N);
if strcmp(naca,'4digit') % 4-digit equation
    m = varargin{1};
    
    if m ~= 0
        yc(1:I) = (m/p^2)*(2*p.*x(1:I) - x(1:I).^2);
        yc(I+1:end) = (m/(1-p)^2)*(1 - 2*p + 2*p.*x(I+1:end) - x(I+1:end).^2);
        
        dyc_dx(1:I) = (2*m/p^2)*(p-x(1:I));
        dyc_dx(I+1:end) = (2*m/(1-p)^2)*(p-x(I+1:end));
    end
    
elseif strcmp(naca,'5digit')  % 5-digit equation 
    q = varargin{1};
    cl = varargin{2};
    
    if q == 0 % Not reflex
        [r,k1] = constants(p,q,cl);
        yc(1:I) = (k1/6)*(x(1:I).^3 - 3*r.*x(1:I).^2 + x(1:I).*r^2*(3-r));
        yc(I+1:end) = (k1*r.^3)./6 .* (1 - x(I+1:end));
        
        dyc_dx(1:I) = (k1/6).* (3*x(1:I).^2 - 6*r.*x(1:I) + r^2 *(3-r));
        dyc_dx(I+1:end) = -k1*r.^3 ./ 6;
        
    elseif q == 1 % reflected camber
        [r,k1,k2] = constants(p,q,cl);
        yc(1:I) = (k1/6)*((x(1:I)-r ).^3 - x(1:I)*k2*(1-r)^3 - ...
            x(1:I)*r^3 + r^3);
        yc(I+1:end) = (k1/6)*(k2*(x(I+1:end) - r).^3 - ...
            x(I+1:end).*k2*(1-r)^3 - x(I+1:end).*r^3 + r^3);
        
        dyc_dx(1:I) = (k1/6).* (3.*(x(1:I)-r).^2 - k2*(1-r)^3 - r^3);
        dyc_dx(I+1:end) = (k1/6).*(3.*k2.*(x(I+1:end) - r).^2 - k2.*(1 - r)^3 - r^3);
    else
        error('not a valid NACA 5 digit airfoil');
    end
    
    
end

if dyc_dx ~= zeros(1,N);
    zeta = atan(dyc_dx);
else
    zeta = 0;
end
end

%% get the constanst for different airfoil
function [m,k1,k2] = constants(p,q,cl)

if q == 0 && nargout == 2
    A = [.05 0.0580 361.400 ;
         .10 .1260  51.640 ;
         .15 .2025  15.957 ;
         .20 .2900  6.643  ;
         .25 .3910  3.230] ;
     
     I = find(p == A(:,1),1);

    if isempty(I)
        i = round(I/5);
        A1 = A(i,:);

        m = (p/A1(1))*A1(2);
        k1 = (p/A1(1))*A1(3);
    else
        m = A(I,2);
        k1 = A(I,3);
    end
elseif q == 1 && nargout == 3
    A = [.10 .1300 51.990 .000764;
         .15 .2170 15.793 .00677;
         .20 .3180 6.520  .0303;
         .25 .4410 3.191  .1355];
     
     I = find(p == A(:,1),1);

    if isempty(I)
        i = round(I/5);
        A1 = A(i,:);

        m = (p/A1(1))*A1(2);
        k1 = (p/A1(1))*A1(3);
    else
        m = A(I,2);
        k1 = A(I,3);
        k2 = A(I,4);
    end
end

% scaling 
if cl ~= .3 && nargout == 2
    m = m*cl/.3;
    k1 = k1*cl/.3;
elseif cl ~= .3 && nargout == 3
    m = m*cl/.3;
    k1 = k1*cl/.3;
    k2 = k2*cl/.3;
end

end