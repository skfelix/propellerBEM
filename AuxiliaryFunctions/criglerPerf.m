function analysis = criglerPerf(inputs,geom,calculate_polars,range,beta_ref_deg,x_beta_ref,varargin)

if nargin < 7
    regime = 'dynamic';
else
    regime = varargin{1};
end

%% Inputs 
V = range;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
unit_system = inputs.unit_system;
x_hub = geom.x_hub;

% Original geometry
x_original = geom.x; % x = r/R
c_original  = geom.c;
beta_original = geom.beta;
betadeg_original = geom.betadeg;

% Transpose arrays if required
[row,~] = size(x_original);
if row > 1
    x_original          = x_original.';
    c_original          = c_original.';
    beta_original       = beta_original.';
    betadeg_original    = betadeg_original.';
end

% Interpolate geom to default x = r/R vector
x = [0.2:0.05:0.95, 0.975, 0.99];
c       = interp1(x_original,c_original,x,'linear','extrap');
beta_interp    = interp1(x_original,beta_original,x,'linear','extrap');
betadeg_interp = interp1(x_original,betadeg_original,x,'linear','extrap');

% Matching beta ref for r = 0.75R or 0.70R - changes blade pitch
delta_betadeg = beta_ref_deg - interp1(x,betadeg_interp,x_beta_ref);
delta_beta = deg2rad(delta_betadeg);
beta = beta_interp + delta_beta;
betadeg = betadeg_interp + delta_betadeg;

% Calculated vars
R = D/2;
r = x*R;
r_hub = x_hub*R;
sigma = B*c./(2*pi*r);
AF = 1e5/(32*R^5)*trapz(x*R, c.*(x*R).^3);
TAF = B*AF;
sigma75 = B*interp1(x,c,0.75)/(2*pi*0.75*R);
rpm = 60*n;
omega = n*2*pi;
S = pi*R^2; % Projected area of the helix

% Atmosphere
T_offset = 0;
[rho,~,~,asound, mu] = ISA(h,T_offset,unit_system);

% Setup 
nP = 'singleRotation'; % 'singleRotation' or 'dualRotation' for contrarrotative props
interpMode = 'linear'; % 'linear' or 'spline


%% Polars

Vref = mean(V);
W0 = sqrt(Vref.^2 + (omega.*x.*R).^2);
W075 = interp1(x,W0,0.75);
Re = rho*W0.*c/mu;
Re75 = interp1(x,Re,0.75);
Mach_ref = W075/asound;

% Limit Mach in 0.45 for xfoil 
if Mach_ref > 0.45
    Mach_ref = 0.45;
end

% Calculate 2D blade section airfoil aerodynamics with xfoil
if (calculate_polars == 1)
    fprintf('======= Running XFoil =======\n')
    [pol,foil] = xfoil(geom.airfoil,-4,20,0.25,Re75,Mach_ref,'oper iter 150','ppar n 200','ppar t 1');
    save('polar.mat','pol','foil');
    fprintf('======= XFoil Finished =======\n')
elseif calculate_polars == 0
    fprintf('======= Polars = 0 =======\n')
    load('polar.mat','pol','foil');
    fprintf('======= XFoil Loaded =======\n')
end

%% Main loop

% Initiate arrays
Cq  = zeros(1,length(V));
Cp  = zeros(1,length(V));
Ct  = zeros(1,length(V));
eta = zeros(1,length(V));
Jw  = zeros(1,length(V));
CL  = zeros(1,length(x));
CD  = zeros(1,length(x));

% Initial guess for interference velocity
wbar = 0.1;

fprintf('%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n','J','CT','CQ','CP','eta','iter')
for i = length(V):-1:1 % Start from smaller speed 

    J(i) = V(i)/(n*D);
    phi0 = atan(V(i)./(omega*r));
    phi = atan(J(i)*(1+0.5*wbar)./(pi*x));
    phi_deg = rad2deg(phi);
    W0 = sqrt(V(i)^2 + (omega.*x.*R).^2);

    % Loop for resultant blade angle (phi) convergence
    error = 10e6;
    count = 1;
    while error > 1e-3

        alpha = beta - phi;
        alphadeg = rad2deg(alpha);

        for j = 1:length(r)
            % Select polar extrapolation method
            pol360 = spera360(pol,foil,geom,1);
            % pol360 = postStallDuSelig(pol,R,V(i),n,r(j),c(j));
            % pol360 = postStallCorriganSchillings(pol,r(j),c(j));
            % pol360 = viternaCorrigan360(pol3D,foil);
            CL(j) = interp1(pol360.alpha,pol360.CL,alphadeg(j),'linear','extrap');
            CD(j) = interp1(pol360.alpha,pol360.CD,alphadeg(j),'linear','extrap');
        end
        
        % Aerodynamic coeffs Reynolds and Mach correction
        %%% CD Mach Reynolds Correction, using Re75 as reference.
        W = V(i)*(1+0.5*wbar*cos(phi).^2)./sin(phi);
        Re = rho*W.*c/mu;
        Mach = W/asound;
        Mach(Mach>0.6) = 0.6; % Limit mach

        % CD = CD.*(Re75./Re).^0.20; % Relation from turbulent BL eq: Cf = 0.074/Re^(1/5) Hernandez & Crespo (1987)
        CD = CD.*(Re75./Re).^0.11; % Relation from experimental data in NACA TR 586. More reliable because BL is not always turbulent on airfoil
        CD = CD.*sqrt(1-Mach_ref.^2)./sqrt(1-Mach.^2);
        CL = CL.*sqrt(1-Mach_ref.^2)./sqrt(1-Mach.^2);
        
        % Decomposed aero coeffs in local blade section axis
        C_axial      = CL.*cos(phi) - CD.*sin(phi);
        C_tangential = CL.*sin(phi) + CD.*cos(phi);
         
        % Differential torque
        dCqdx = 1/8*pi*x.^2*J(i)^2.*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigma.*C_tangential; % From Borst (1973)
        Cq(i) = trapz(x,dCqdx);
        Cp(i) = 2*pi*Cq(i);

        Pct = 64/pi*(n*R/V(i))^3*Cp(i); % Cp to Pc
        options = optimset('Display','off');
        Pc = @(wbar) 2*getKappa(nP, interpMode, B, J(i)*(1+wbar))*wbar*(1+wbar)*(1+get_eKappa(nP, interpMode, B, J(i)*(1+wbar))*wbar);
        wbarnew = fsolve(@(wbar)Pc(wbar)-Pct,wbar,options);
        Jw(i) = J(i)*(1+wbarnew);

        % Calculate new resultant blade angle
        phinew = atan(J(i)*(1+0.5*wbarnew)./(pi*x));

        error = max(abs(wbar-wbarnew));
        
        % Update values - with relaxation factor
        k = 0.5;
        wbar = (1-k)*wbar + k*wbarnew;
        phi = (1-k)*phi + k*phinew;
        
        % Limit loop in 50 iterations
        count = count + 1;
        if count > 50
            break;
        end
    end

    % Calculate propulsive coeffs
    dCtdx = 1/4*pi*x*J(i)^2.*((1+0.5*wbar*cos(phi).^2)./sin(phi)).^2.*sigma.*C_axial;  % From Borst (1973)
    Ct(i) = trapz(x,dCtdx);

    eta(i) = J(i)*Ct(i)/Cp(i);
    T(i) = Ct(i)*(rho*n^2*D^4);
    Q(i) = Cq(i)*(rho*n^2*D^5);
    P(i) = Cp(i)*(rho*n^3*D^5);

    % Print results
    fprintf('%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n',J(i),Ct(i),Cq(i),Cp(i),eta(i),count)
end

% Assemble output struct
analysis.rpm         = rpm;
analysis.n           = n;
analysis.V           = V(Ct~=0);
analysis.J           = J(Ct~=0);
analysis.Ct          = Ct(Ct~=0);
analysis.Cp          = Cp(Ct~=0);
analysis.Cq          = Cq(Ct~=0);
analysis.eta         = eta(Ct~=0);
analysis.T           = T(Ct~=0);
analysis.P           = P(Ct~=0);
analysis.Q           = Q(Ct~=0);
analysis.unit_system = unit_system;
analysis.Re75        = Re75;
analysis.Mach75      = Mach_ref;