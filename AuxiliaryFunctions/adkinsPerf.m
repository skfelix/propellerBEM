function analysis = adkinsPerf(inputs,geom,calculate_polars,range,beta_ref_deg,x_beta_ref)







%% Inputs 
V = range;
n = inputs.n;
D = inputs.D;
B = inputs.B;
h = inputs.h;
unit_system = inputs.unit_system;
x_hub = geom.x_hub;

% Original geometry
x_original = geom.x;
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

%% Polars

Vref = mean(V);
W0 = sqrt(Vref.^2 + (omega.*x.*R).^2);
W075 = interp1(x,W0,0.75);
Re = rho*W0.*c/mu;
Re75 = interp1(x,Re,0.75);
Mach_ref = W075/asound;
% Mach_ref = 0;

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
J   = zeros(1,length(V));
Cp  = zeros(1,length(V));
Cq  = zeros(1,length(V));
Ct  = zeros(1,length(V));
eta = zeros(1,length(V));
CL  = zeros(1,length(x));
CD  = zeros(1,length(x));


fprintf('%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n','J','CT','CQ','CP','eta','iter')
for i = length(V):-1:1

    J(i) = V(i)/(n*D);
    omega = n*2*pi;
    lambda = V(i)/(omega*R);
    chi = omega.*r/V(i);
    
    phi0 = atan(V(i)./(omega*r));    
    W0 = sqrt(V(i)^2 + (omega.*x.*R).^2);
    
    % Initial Guesses
    W = W0; 
    phi = phi0;
    
    % Loop for resultant blade angle (phi) convergence
    error = 10e6;
    count = 1;
    while error > 1e-4
        
        phit = atan(tan(phi(end)).*x(end));
        f  = 0.5*B*(1-x)./sin(phit); % adkins
        F  = 2/pi*acos(exp(-f));
        fh = 0.5*B*(r-r_hub)./(r_hub*sin(phi));
        Fh = 2/pi*acos(exp(-fh));
        F  = F.*Fh;
        
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
        Re = rho*W.*c/mu;
        Mach = W/asound;
        Mach(Mach>0.60) = 0.60; % Limit mach

        % CD = CD.*(Re75./Re).^0.20; % Relation from turbulent BL eq: Cf = 0.074/Re^(1/5) Hernandez & Crespo (1987)
        CD = CD.*(Re75./Re).^0.11; % Relation from experimental data in NACA TR 586. More reliable because BL is not always turbulent on airfoil
        CD = CD.*sqrt(1-Mach_ref.^2)./sqrt(1-Mach.^2);
        CL = CL.*sqrt(1-Mach_ref.^2)./sqrt(1-Mach.^2);        
        
        % Decomposed aero coeffs in local blade section axis
        C_axial      = CL.*cos(phi) - CD.*sin(phi); 
        C_tangential = CL.*sin(phi) + CD.*cos(phi);
        
        % Interference velocities
        Ka = C_axial./(4*sin(phi).^2);
        Kb = C_tangential./(4*cos(phi).*sin(phi));
        a = sigma.*Ka./(F - sigma.*Ka);
        b = sigma.*Kb./(F + sigma.*Kb);
        
        % Set tip phi value as zero (as per adkins)
        if phi(end) < 0
            phi(end) = deg2rad(1);
        end
        
        % Calculate new resultant blade angle
        phinew = atan(V(i).*(1+a)./(omega.*r.*(1-b)));
        
        error = max(abs(phi-phinew));
        
        % Update values - with relaxation factor
        k = 0.5; % Relaxation factor
        phi = (1-k)*phi + k*phinew;

        count = count + 1;
        % Limit loop in 50 iterations
        if count > 50
            break;
        end
    end
    
    % Calculate propulsive coeffs
    W = V(i).*(1+a)./sin(phi);
    
    dTdx = 0.5*rho.*W.^2*B.*c.*C_axial;
    T(i) = trapz(r,dTdx);
    dQdx = 0.5*rho.*W.^2*B.*c.*C_tangential.*r;
    Q(i) = trapz(r,dQdx);
    P(i) = omega*Q(i); % P = 2*pi*n*Q = omega*Q
    Ct(i) = T(i)./(rho*n^2*D^4);
    Cq(i) = Q(i)./(rho*n^2*D^5);
    Cp(i) = P(i)./(rho*n^3*D^5);
    eta(i) = J(i)*Ct(i)/Cp(i);

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