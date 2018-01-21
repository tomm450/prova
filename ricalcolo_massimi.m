for i = 1:9
    
    fprintf('i = %d \n',i)
    x_t{i} = linspace(x(1)-1,x(1)+1,4);
    
    
end
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

