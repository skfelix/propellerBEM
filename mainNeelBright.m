%% MAIN PROP - ALL BETA
clc; clear; fclose all; 

% Add folders to path  
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%%%% Neel & Bright (1950) - imperial units
load('neelBright_data.mat')
propData = neelBright;
Wdot = propData.perf.beta21.Wdot;
n = propData.perf.beta21.n;
D = propData.geom.D;
B = propData.geom.B;
h = 0;
unit_system = propData.geom.unit_system;
propData.geom.airfoil = 'clark-y.dat';

inputs.Wdot = Wdot; inputs.n = n; inputs.D = D; inputs.B = B; inputs.h = h; inputs.unit_system = unit_system;

temp_offset = 0; % 10 for ISA+10, or -10 for ISA-10 for example
rho = ISA(h,0,unit_system);
polars  = 0;
x        = propData.geom.x;
betadeg    = propData.geom.betadeg;

Cp1 = Wdot/(rho*n^3*D^5);

run = 0;
% Performance
if run 
    % BETA75 = 20 deg
    V = 140:10:230;
    beta_ref_deg = propData.perf.beta21.betadeg;
    x_beta_ref = propData.perf.beta21.x_beta;
    an1   = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an2   = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 25 deg
    V = 160:10:280;
    beta_ref_deg = propData.perf.beta25.betadeg;
    x_beta_ref = propData.perf.beta25.x_beta;
    an3   = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an4   = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 30 deg
    V = 200:10:360;
    beta_ref_deg = propData.perf.beta31.betadeg;
    x_beta_ref = propData.perf.beta31.x_beta;
    an5   = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an6   = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    save('results/neelBright_results.mat','an1','an2','an3','an4','an5','an6');
else
    load('results/neelBright_results.mat','an1','an2','an3','an4','an5','an6');
end

%% Plots

nfig = 1;
lgd= {
    '$\beta_{75}=21^\circ$ - Crigler'
    '$\beta_{75}=25^\circ$ - Crigler'
    '$\beta_{75}=31^\circ$ - Crigler'
    '$\beta_{75}=21^\circ$ - Adkins'
    '$\beta_{75}=25^\circ$ - Adkins'
    '$\beta_{75}=31^\circ$ - Adkins'
    'Experimental'
    };
title = strcat('Neel & Bright (1950)');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 


% Perf
% nfig = nfig+1;
plotVelocity = 0;
figure(nfig);
perf1 = {an1,an3,an5};
perf2 = {an2,an4,an6};
group = 1;
plotPerf(perf1,title,lgd,group,plotVelocity)
group = 2;
plotPerf(perf2,title,lgd,group,plotVelocity)
hold on
subplot(211)
plot(propData.perf.beta21.J_Ct,propData.perf.beta21.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta25.J_Ct,propData.perf.beta25.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta31.J_Ct,propData.perf.beta31.Ct,'ok','MarkerFaceColor','k')
subplot(212)
plot(propData.perf.beta21.J_Cp,propData.perf.beta21.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta25.J_Cp,propData.perf.beta25.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta31.J_Cp,propData.perf.beta31.Cp,'ok','MarkerFaceColor','k')