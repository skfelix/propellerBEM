function pol_extended = viternaCorrigan360(pol,foil)
% [pol,foil] = xfoil('clark-y.dat',-10:0.5:20,2.3e6,0,'oper iter 150');
pol_extended = pol;

t_c = foil.max_thickness;
rLE = foil.rLE;

% Subscript s: at stall
CL_s = max(pol.CL);
alpha_s = interp1(pol.CL,pol.alpha,CL_s);
CD_s = interp1(pol.alpha,pol.CD,alpha_s);

% Constants
CD90 = 2.0772-3.978*rLE;% $ CDmax
B1 = CD90;
B2 = CD_s-CD90*sind(alpha_s).^2./cosd(alpha_s);
A1 = B1/2;
A2 = (CL_s-CD90*sind(alpha_s)*cosd(alpha_s))*sind(alpha_s)./cosd(alpha_s).^2;

alpha = -10:0.5:120;

for i = 1:length(alpha)
    if alpha(i) <= max(pol.alpha)
        CL(i) = interp1(pol.alpha,pol.CL,alpha(i),'linear','extrap');
        CD(i) = interp1(pol.alpha,pol.CD,alpha(i),'linear','extrap');
        if isnan(CL(i))
            aaa= 1;
        end
    elseif alpha(i) > max(pol.alpha)               
        CL(i) = A1*sind(2*alpha(i)) + A2*cosd(alpha(i)).^2./sind(alpha(i));
        CD(i) = B1*sind(alpha(i)).^2 + B2*cosd(alpha(i));
    end
end


pol_extended.alpha = alpha;
pol_extended.CL = CL;
pol_extended.CD = CD;

% figure(1)
% subplot(121)
% hold on
% plot(pol.alpha,pol.CL)
% plot(alpha,CL,'--')
% 
% subplot(122)
% hold on
% plot(pol.alpha,pol.CD)
% plot(alpha,CD,'--')
% 
% abd = 5;







