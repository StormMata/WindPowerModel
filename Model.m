clearvars -except AeroDyn A Output Hannah HannahSS SSs SSsn speed uRel vRel Hphi aoa Hchord Htwist Hpower HW sweep
clc

speed = -0.1:0.01:0.7;

% for s = 1:length(speed)

inc     = 0;
Palpha  = 0.0;
% Palpha  = speed(s);

z       = (104+63:-1:104-63)';
% beta    = linspace(12.6,-12.6,length(z))';
beta    = flip(0 + inc * (0:length(z)-1)');
H2      = ceil(length(z)/2);
B2      = beta(H2);
% B1      = beta(T2 + 1)
% mid     = abs(B1 + B2)/2
beta    = beta - B2;

yaw     = 0.0;
azimuth = linspace(0,2 * pi,360);
% U       = 11.4;% * (1 - 1/3);
r       = AeroDyn.r + 1.5;
Omega   = (12.1 * 2 * pi)/(60);
rho     = 1.225;
A_po    = 0;                                                            % Pitch angle                   [rad]
B       = 3;                                                            % Number of blades              [-]
R_hub   = 1.5;                                                          % Hub radius                    [m]
R_rot   = 63;
% a       = Output(:,2);
% a       = [a' a(end) a(end)];
a       = 1/3 * ones(1, length(r));

for j = 1:length(azimuth)

    for i = 1:length(r)

%         U      = 11.4 * ((104 + r(i) * cos(azimuth(j)))/(104))^0.2;

        zloc   = round(104 + r(i) * cos(azimuth(j)),0);

        betaind= find(zloc==z);

        betar  = deg2rad(beta(betaind));

        U      = 11.4 * ((104 + r(i) * cos(azimuth(j)))/(104))^Palpha;

        U_axi(i)= (U*(1-a(i)) * cos(betar * sin(azimuth(j))) * cos(betar * cos(azimuth(j))));

        U_tan(i)= (Omega * r(i) - U*(1-a(i)) * cos(betar * sin(azimuth(j))) * sin(betar * cos(azimuth(j))));
    
        Wm(i)  = (U_axi(i)/11.4)^2 + (U_tan(i)/11.4)^2;

        W(i)   = U_axi(i)^2 + U_tan(i)^2;
    
        Phi(i) = atan2(U_axi(i), U_tan(i));
    
        c      = AeroDyn.Chord(i);
    
        alpha(i)  = rad2deg(Phi(i)) - AeroDyn.Twist(i);

        CL(i) = interp1(A.(sprintf('%s',AeroDyn.AeroID(i))).Alpha,A.(sprintf('%s',AeroDyn.AeroID(i))).CL,alpha(i),'linear');
        CD(i) = interp1(A.(sprintf('%s',AeroDyn.AeroID(i))).Alpha,A.(sprintf('%s',AeroDyn.AeroID(i))).CD,alpha(i),'linear');

        dQ(i)  = 3 * 0.5 * rho * c * W(i) * r(i) * (CL(i) * sin(Phi(i)) - CD(i) * cos(Phi(i)));
    
        dQm(i) = c/63 * r(i)/63 * Wm(i) * (CL(i) * sin(Phi(i)) - CD(i) * cos(Phi(i)))
    end

    dQ    = trapz(r,dQ);

    dP(j) = Omega * dQ;

end

trapz(azimuth,dP)/(2*pi)
% sweep(s) = trapz(azimuth,dP)/(2*pi);

% end

% figure;
% subplot(2,2,[1 2])
%     plot(r,U_axi,'b','LineWidth',1.5)
%         hold on
%     plot(r,uRel,'r','LineWidth',1.5)
%         title('Axial Velocity')
%         legend('Storm','Hannah','Location','northeast')
%         ylim([7.59999998 7.60000007])
%         xlabel('r (m)')
%         ylabel('uRel (m/s)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (U_axi - uRel)./((U_axi + uRel)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
%         ylim([-6e-7 -4e-7])
% 
% figure;
% subplot(2,2,[1 2])
%     plot(r,U_tan,'b','LineWidth',1.5)
%         hold on
%     plot(r,vRel,'r','LineWidth',1.5)
%         title('Tangential Velocity')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('r (m)')
%         ylabel('vRel (m/s)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (U_tan - vRel)./((U_tan + vRel)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
% 
% figure;
% subplot(2,2,[1 2])
%     plot(r,Phi,'b','LineWidth',1.5)
%         hold on
%     plot(r,Hphi,'r','LineWidth',1.5)
%         title('Inflow Angle (\phi) ')
%         legend('Storm','Hannah','Location','northeast')
%         xlabel('r (m)')
%         ylabel('\phi (rad)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (Phi - Hphi)./((Phi + Hphi)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
% 
% figure;
% subplot(2,2,[1 2])
%     plot(r,alpha,'b','LineWidth',1.5)
%         hold on
%     plot(r,aoa,'r','LineWidth',1.5)
%         title('Angle of Attack (\alpha) ')
%         legend('Storm','Hannah','Location','northeast')
%         xlabel('r (m)')
%         ylabel('\alpha (deg)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (alpha - aoa)./((alpha + aoa)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
% 
% figure;
% subplot(2,2,[1 2])
%     plot(r,AeroDyn.Chord,'b','LineWidth',1.5)
%         hold on
%     plot(r,Hchord,'r','LineWidth',1.5)
%         title('Chord')
%         legend('Storm','Hannah','Location','northeast')
%         xlabel('r (m)')
%         ylabel('Chord (m)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (AeroDyn.Chord - Hchord')./((AeroDyn.Chord + Hchord')/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
% % 
% figure;
% subplot(2,2,[1 2])
%     plot(r,AeroDyn.Twist,'b','LineWidth',1.5)
%         hold on
%     plot(r,Htwist,'r','LineWidth',1.5)
%         title('Twist')
%         legend('Storm','Hannah','Location','northeast')
%         xlabel('r (m)')
%         ylabel('Twist (deg)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (AeroDyn.Twist - Htwist')./((AeroDyn.Twist + Htwist')/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])

% figure;
% subplot(2,2,[1 2])
%     plot(r,W,'b','LineWidth',1.5)
%         hold on
%     plot(r,HW,'r','LineWidth',1.5)
%         title('U_{rel}^2')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('r (m)')
%         ylabel('U_{rel}^2 (m^2/s^2)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (W - HW)./((W + HW)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
% 
figure;
subplot(2,2,[1 2])
    plot(speed,sweep,'b','LineWidth',1.5)
        hold on
    plot(speed,Hpower,'r','LineWidth',1.5)
        title('Power')
        legend('Storm','Hannah','Location','northwest')
        xlabel('Speed Shear Exponent (\alpha)')
        ylabel('Power (MW)')
        xlim([min(speed) max(speed)])
subplot(2,2,[3 4])
    plot(speed, 100 * (sweep - Hpower)./((sweep + Hpower)/2),'LineWidth',1.5)
        xlabel('Speed Shear Exponent (\alpha)')
        ylabel('Percent Difference (%)')
        xlim([min(speed) max(speed)])
