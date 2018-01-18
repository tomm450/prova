clc
clear all
close all

%% CL-alpha
% report 824, pagina 445 (188 nel pdf)
% immagine a base green 
% libero due canali
% % decommentare se si deve risvuotare canali R e B 
% in_1 = imread('./1.png');
% in_2 = imread('./2.png');
% 
% in_1(:,:,1) = zeros(size(in_1(:,:,1)));
% in_2(:,:,1) = zeros(size(in_2(:,:,1)));
% in_1(:,:,3) = zeros(size(in_1(:,:,3)));
% in_2(:,:,3) = zeros(size(in_2(:,:,3)));
% 
% imwrite(in_1,'1green.png');
% imwrite(in_2,'2green.png');

% blu per assi [0 0 80] per i 5 punti (origine e cardinali)
% rosso per curve (Cl-alpha, Cm-alpha, Cd-Cl)
% - Re = 3M rappresentato con 'o', rgb = [ 80 0 0]
% - Re = 6M rappresentato con 's', rgb = [160 0 0]
% - Re = 9M rappresentato con 'd', rgb = [240 0 0]
% NB: punti ricalcati a "mano" con pinta

% Dati prova %% PRESSURE TUNNEL, rho non nota
% Chord = 24*2.54/100; % m
% Nu    = 0.000018375; %kg/m s
% Rho   = 1.225;       %kg/m3
% sos   = 340;         % m/s
% 
% RE   = [3; 6; 9]*1000000;
% 
% U_oo = RE*Nu./(Rho*Chord);
% M_oo = U_oo./sos;


nptNew = 100;
plt_ord_raw = {'ro--','cs--','bd--'};
plt_ord     = {'r-','c-','b-'};

%% assiANDcl & cm.png
% origine
o = [0 0 80];
% asse x 
x_fs = [-32,32];     DA_FS = x_fs(2)-x_fs(1);
x_mk = [0 0 240];
% asse y
cl_fs = [-2,3.6];    DCL_FS = cl_fs(2)-cl_fs(1);
cm_fs = [-0.5,0.9];  DCM_FS = cm_fs(2)-cm_fs(1);

[x,y] = lettura_canali2('assiANDcl.png',o);
cross = [x;y];

[left,ileft]   = min(cross(1,:));
[right,iright] = max(cross(1,:));
[top,itop]     = min(cross(2,:));
[bot,ibot]     = max(cross(2,:));

[or,ior ] = min(sum( [ileft iright itop ibot]' == [1 2 3 4 5]));

LR = cross(:,iright) - cross(:,ileft);
TB = cross(:,ibot) - cross(:,itop);

% inclinazione asse
a = mean([atan2d(LR(2),LR(1)) , atan2d(-TB(1),TB(2))]);
a = -a;
Rmat = [ cosd(a) -sind(a) ; ...
         sind(a)  cosd(a)];
     
% dimensioni px
Xpx = norm(LR);
Ypx = norm(TB);

%% Cl-Alpha
figure(1);
cl_a_raw = cell(1,3);
cl_a     = cell(1,3);

for i = 1:3
   
   [x_temp,y_temp] = lettura_canali2('assiANDcl.png',[80*i 0 0]);
   % creo matrice
   cross_temp = [x_temp;y_temp];
   
   % ordino secondo x
   [~,idx] = sort(cross_temp(1,:));
   cross_temp = cross_temp(:,idx);
   
   % tolgo offset
   cross_temp = cross_temp - cross(:,ior);
   
   % ruoto per correggere orientamento
   cross_temp = Rmat*cross_temp;
   
   % ribalto asse
   cross_temp(1,:) =  cross_temp(1,:)*(DA_FS/Xpx);
   cross_temp(2,:) = -cross_temp(2,:)*(DCL_FS/Ypx);
   
   a_fine = linspace(min(cross_temp(1,:)),max(cross_temp(1,:)),nptNew);
   c_fine = interp1(cross_temp(1,:),cross_temp(2,:),a_fine);
   
   cl_a_raw{i} = cross_temp; 
   cl_a{i}     = [a_fine;c_fine]; 
   
   plot(cross_temp(1,:),cross_temp(2,:),plt_ord_raw{i})
   hold on
   plot(a_fine,c_fine,plt_ord{i})
   grid on
   title('red = Re 3M, cian = Re 6M, blue = Re 9M')
end

%% Cm-Alpha

cm_a_raw = cell(1,3);
cm_a     = cell(1,3);
figure(2)
for i = 1:3
   
   [x_temp,y_temp] = lettura_canali2('cm.png',[80*i 0 0]);
   % creo matrice
   cross_temp = [x_temp;y_temp];
   
   % ordino secondo x
   [~,idx] = sort(cross_temp(1,:));
   cross_temp = cross_temp(:,idx);
   
   % tolgo offset
   cross_temp = cross_temp - cross(:,ior);
   
   % ruoto per correggere orientamento
   cross_temp = Rmat*cross_temp;
   
   % ribalto asse
   cross_temp(1,:) =  cross_temp(1,:)*(DA_FS/Xpx);
   cross_temp(2,:) = -cross_temp(2,:)*(DCM_FS/Ypx);
   
   a_fine = linspace(min(cross_temp(1,:)),max(cross_temp(1,:)),nptNew);
   c_fine = interp1(cross_temp(1,:),cross_temp(2,:),a_fine);
   
   cm_a_raw{i} = cross_temp; 
   cm_a{i}     = [a_fine;c_fine]; 
   
   plot(cross_temp(1,:),cross_temp(2,:),plt_ord_raw{i})
   hold on
   grid on
   plot(a_fine,c_fine,plt_ord{i})
   grid on
   title('red = Re 3M, cian = Re 6M, blue = Re 9M')
end


%% assiANDpol 
% origine
o = [0 0 80];
% asse x 
cl_fs = [-1.6,1.6];    DCL_FS = cl_fs(2)-cl_fs(1);
% asse y
cd_fs = [-0.02,0.024]; DCD_FS = cd_fs(2)-cd_fs(1);

[x,y] = lettura_canali2('assiANDpol.png',o);
cross = [x;y];

[left,ileft]   = min(cross(1,:));
[right,iright] = max(cross(1,:));
[top,itop]     = min(cross(2,:));
[bot,ibot]     = max(cross(2,:));

[or,ior ] = min(sum( [ileft iright itop ibot]' == [1 2 3 4 5]));

LR = cross(:,iright) - cross(:,ileft);
TB = cross(:,ibot) - cross(:,itop);

% inclinazione asse
a = mean([atan2d(LR(2),LR(1)) , atan2d(-TB(1),TB(2))]);
a = -a;
Rmat = [ cosd(a) -sind(a) ; ...
         sind(a)  cosd(a)];
     
% dimensioni px
Xpx = norm(LR);
Ypx = norm(TB);

% Cl-Cd
figure(3);
cl_cd_raw = cell(1,3);
cl_cd     = cell(1,3);

for i = 1:3
   
   [x_temp,y_temp] = lettura_canali2('assiANDpol.png',[80*i 0 0]);
   % creo matrice
   cross_temp = [x_temp;y_temp];
   
   % ordino secondo x
   [~,idx] = sort(cross_temp(1,:));
   cross_temp = cross_temp(:,idx);
   
   % tolgo offset
   cross_temp = cross_temp - cross(:,ior);
   
   % ruoto per correggere orientamento
   cross_temp = Rmat*cross_temp;
   
   % ribalto asse
   cross_temp(1,:) =  cross_temp(1,:)*(DCL_FS/Xpx);
   cross_temp(2,:) = -cross_temp(2,:)*(DCD_FS/Ypx);
   
   a_fine = linspace(min(cross_temp(1,:)),max(cross_temp(1,:)),nptNew);
   c_fine = interp1(cross_temp(1,:),cross_temp(2,:),a_fine);
   
   cl_cd_raw{i} = cross_temp; 
   cl_cd{i}     = [a_fine;c_fine]; 
   
   plot(cross_temp(1,:),cross_temp(2,:),plt_ord_raw{i})
   hold on
   grid on
   plot(a_fine,c_fine,plt_ord{i})
   grid on
   title('red = Re 3M, cian = Re 6M, blue = Re 9M')
   
end



%% subroutine 
function [x,y] = lettura_canali2(nome_img,c2r)

%c2r = canale da leggere

x = [];
y = [];

a = imread(nome_img);

if c2r(2) == 0 && c2r(3) == 0
    c = 1;
elseif c2r(1) == 0 && c2r(3) == 0
    c = 2;
elseif c2r(1) == 0 && c2r(2) == 0
    c = 3;
else 
    error('c2r == [0 0 0] ?');
end
    
nfind = sum(sum(a(:,:,c) == c2r(c))) ;
success = 0;

while success < nfind
    
    % sum(a(:,:,c) == c2r(c))) vettore colonna 
    [p_col, icol ] = max(sum(a(:,:,c) == c2r(c)));
    
    for k = 1:p_col
        [toBeOne, irow ] = max(a(:,icol,c) == c2r(c));
        if toBeOne == 1
            % ok
        x = [x,icol];
        y = [y,irow];
        
        success = success +1;
        a(irow,icol,c) = 0;
        else
            error('toBeOne should be one duh!')
        end
    end
    
    
end


end


