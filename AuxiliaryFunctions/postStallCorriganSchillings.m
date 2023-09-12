function pol3D = postStallCorriganSchillings(pol,r,c)

pol3D = pol;

alpha = pol.alpha;
CL2D = pol.CL;
alpha0 = interp1(CL2D,alpha,0);
CL2D_max = max(CL2D);
alpha_max = interp1(CL2D,alpha,CL2D_max);
linReg = polyfit(alpha(alpha>-5 & alpha<5),CL2D(alpha>-5 & alpha<5),1);         
CLa = linReg(1); % dCL/dalpha 

n = 1;
c_r = c/r;
K = (0.1517/(c_r))^(1/1.084);
delta_alpha = ((K*(c_r)/0.136)^n-1)*(alpha_max-alpha0);

alpha3D = 0*alpha;
CL3D = 0*CL2D;

for i = 1:length(alpha)
    if alpha(i) >= 5
        alpha3D(i) = alpha(i) + delta_alpha;
        CL3D(i) = CL2D(i) + CLa*delta_alpha;
    else
        alpha3D(i) = alpha(i);
        CL3D(i) = CL2D(i);
    end
end

pol3D.alpha = alpha3D;
pol3D.CL = CL3D;

% figure(1)
% hold on
% plot(alpha,CL2D)
% plot(alpha3D,CL3D)
% % plot(pol.alpha,pol.CL)
% % plot(alpha,CL2D)
% 
% 
% abd = 5;
