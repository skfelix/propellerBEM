%% MAIN PROP - ALL BETA
clc; clear; fclose all; 

% Add folders to path  
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));


%%%% Corson & Maynard (1946) - imperial
load('corsonMaynard_data.mat')
propData = corsonMaynard;
Wdot = propData.perf.beta20.Wdot;
n = propData.perf.beta20.n;
D = propData.geom.D;
B = propData.geom.B;
h = 0;
unit_system = propData.geom.unit_system;
propData.geom.airfoil = 'clark-y.dat';

inputs.Wdot = Wdot; inputs.n = n; inputs.D = D; inputs.B = B; inputs.h = h; inputs.unit_system = unit_system;

rho = ISA(h,0,unit_system);
polars  = 0;
x        = propData.geom.x;
betadeg    = propData.geom.betadeg;

Cp1 = Wdot/(rho*n^3*D^5);

run = 1;
% Performance
if run
    % BETA75 = 20 deg
    V = 160:5:275;
    beta_ref_deg = propData.perf.beta20.betadeg;
    x_beta_ref = propData.perf.beta20.x_beta;
    an1 = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an6 = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 25 deg
    V = 200:5:340;
    beta_ref_deg = propData.perf.beta25.betadeg;
    x_beta_ref = propData.perf.beta25.x_beta;
    an2 = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an7 = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 30 deg
    V = 250:5:405;
    beta_ref_deg = propData.perf.beta30.betadeg;
    x_beta_ref = propData.perf.beta30.x_beta;
    an3 = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an8 = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 35 deg
    V = 350:5:510;
    beta_ref_deg = propData.perf.beta35.betadeg;
    x_beta_ref = propData.perf.beta35.x_beta;
    an4 = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an9 = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    % BETA75 = 40 deg
    V = 470:5:600;
    beta_ref_deg = propData.perf.beta40.betadeg;
    x_beta_ref = propData.perf.beta40.x_beta;
    an5 = criglerPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    an10 = adkinsPerf(inputs,propData.geom,polars,V,beta_ref_deg,x_beta_ref);
    save('results/corsonMaynard_results.mat','an1','an2','an3','an4','an5','an6','an7','an8','an9','an10');
else
    load('results/corsonMaynard_results.mat','an1','an2','an3','an4','an5','an6','an7','an8','an9','an10');
end

%% Plots

nfig = 1;
lgd= {
    '$\beta_{75}=20^\circ$ - Crigler'
    '$\beta_{75}=25^\circ$ - Crigler'
    '$\beta_{75}=30^\circ$ - Crigler'
    '$\beta_{75}=35^\circ$ - Crigler'
    '$\beta_{75}=40^\circ$ - Crigler'
    '$\beta_{75}=20^\circ$ - Adkins'
    '$\beta_{75}=25^\circ$ - Adkins'
    '$\beta_{75}=30^\circ$ - Adkins'
    '$\beta_{75}=35^\circ$ - Adkins'
    '$\beta_{75}=40^\circ$ - Adkins'
    'Experimental'
    };
title = strcat('Corson & Maynard (1946)');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 


% Perf
% nfig = nfig+1;
plotVelocity = 0;
figure(nfig);
perf1 = {an1,an2,an3,an4,an5};
perf2 = {an6,an7,an8,an9,an10};
group = 1;
plotPerf(perf1,title,lgd,group,plotVelocity)
group = 2;
plotPerf(perf2,title,lgd,group,plotVelocity)
hold on
subplot(211)
plot(propData.perf.beta20.J_Ct,propData.perf.beta20.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta25.J_Ct,propData.perf.beta25.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta30.J_Ct,propData.perf.beta30.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta35.J_Ct,propData.perf.beta35.Ct,'ok','MarkerFaceColor','k')
plot(propData.perf.beta40.J_Ct,propData.perf.beta40.Ct,'ok','MarkerFaceColor','k')
subplot(212)
plot(propData.perf.beta20.J_Cp,propData.perf.beta20.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta25.J_Cp,propData.perf.beta25.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta30.J_Cp,propData.perf.beta30.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta35.J_Cp,propData.perf.beta35.Cp,'ok','MarkerFaceColor','k')
plot(propData.perf.beta40.J_Cp,propData.perf.beta40.Cp,'ok','MarkerFaceColor','k')