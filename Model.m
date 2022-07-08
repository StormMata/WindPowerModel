clearvars -except dir hz hb z beta afind diff Storm P Q AeroDyn A Output Hannah HannahSS SSs SSsn speed uRel vRel Hphi aoa Hchord Htwist Hpower HW sweep Hcl Hcd Hr t0 tn2  test
clc

%% Inputs

load('Turbine.mat')                                                         % Load turbine parameters

load('Airfoils.mat')                                                        % Load airfoil parameters

load('WindProfiles.mat')                                                    % Load atmospheric conditions


speed = -0.1:0.1:0.7;

% for s = 1:length(speed)
%     s

inc     = 0.0;
Palpha  = 0.0;
% Palpha  = speed(s);

z       = (T.Hub+63:-1:T.Hub-63)';
beta    = flip(0 + inc * (0:length(z)-1)');
H2      = ceil(length(z)/2);
B2      = beta(H2);
beta    = beta - B2;

yaw     = 0.0;
azimuth = linspace(0,2 * pi,50);
r       = linspace(1.5,63,100);
Omega   = (12.1 * 2 * pi)/(60);
rho     = 1.225;
A_po    = 0;                                                            % Pitch angle                   [rad]
B       = 3;                                                            % Number of blades              [-]
R_hub   = 1.5;                                                          % Hub radius                    [m]
R_rot   = 63;
a       = 1/3 * ones(1, length(r));

for j = 1:length(azimuth)

    for i = 1:length(r)

        zloc           = round(T.Hub + r(i) * cos(azimuth(j)),0);

        betar          = deg2rad(interp1(z,beta,zloc,'linear','extrap'));

%         U_Hub     = interp1(I.Heights,I.SpeedShear(:,i),T.Hub,'linear');

        U(j,i)         = 11.4 * ((T.Hub + r(i) * cos(azimuth(j)))/(T.Hub))^Palpha;

        U_axi(j,i)     = (U(j,i)*(1-a(i)) * cos(betar * sin(azimuth(j))) * cos(betar * cos(azimuth(j))));

        U_tan(j,i)     = (Omega * r(i) - U(j,i)*(1-a(i)) * cos(betar * sin(azimuth(j))) * sin(betar * cos(azimuth(j))));
    
        W(j,i)         = U_axi(j,i)^2 + U_tan(j,i)^2;
    
        Phi(j,i)       = atan2(U_axi(j,i), U_tan(j,i));
    
        c(j,i)         = interp1(AeroDyn.r,AeroDyn.Chord,r(i),'linear','extrap');

        Twist(j,i)     = interp1(AeroDyn.r,AeroDyn.Twist,r(i),'linear','extrap');
    
        alpha(j,i)     = rad2deg(Phi(j,i)) - Twist(j,i);

        AeroIndex(j,i) = interp1(AeroDyn.r,AeroDyn.AeroIndex,r(i),'nearest','extrap');      % Interpolated airfoil index        [-]

        CL(j,i)        = interp1(A.(sprintf('A%i',AeroIndex(j,i))).Alpha,A.(sprintf('A%i',AeroIndex(j,i))).CL,alpha(j,i),'linear');
        CD(j,i)        = interp1(A.(sprintf('A%i',AeroIndex(j,i))).Alpha,A.(sprintf('A%i',AeroIndex(j,i))).CD,alpha(j,i),'linear');

        dQ(j,i)        = 3 * 0.5 * rho * c(j,i) * W(j,i) * r(i) * (CL(j,i) * sin(Phi(j,i)) - CD(j,i) * cos(Phi(j,i)));

    end

end

% end

P = dQ * Omega;

trapz(r,trapz(azimuth,P,1))/(1e6*2*pi)

[Ra Th] = meshgrid(r,azimuth);

X = Ra.*sin(Th);
Y = Ra.*cos(Th);

P = figure;
surf(X,Y+T.Hub,dQ,'LineStyle','none')
xlabel('y (m)')
ylabel('z (m)')
axis equal
view(gca,0,90)


% Storm = [sweep; Storm]

% figure;
% subplot(2,2,[1 2])
%     plot(r,U_axi/11.4,'b','LineWidth',1.5)
%         hold on
%     plot(r,uRel,'r','LineWidth',1.5)
%         title('Axial Velocity')
%         legend('Storm','Hannah','Location','northeast')
%         ylim([0.66666664 0.66666668])
%         xlabel('r (m)')
%         ylabel('uRel (m/s)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (U_axi/11.4 - uRel)./((U_axi/11.4 + uRel)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
%         ylim([-6e-7 -4e-7])
% 
% figure;
% subplot(2,2,[1 2])
%     plot(r,U_tan/11.4,'b','LineWidth',1.5)
%         hold on
%     plot(r,vRel,'r','LineWidth',1.5)
%         title('Tangential Velocity')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('r (m)')
%         ylabel('vRel (m/s)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (U_tan/11.4 - vRel)./((U_tan/11.4 + vRel)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])

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
%     plot(r,c/63,'b','LineWidth',1.5)
%         hold on
%     plot(r,Hchord,'r','LineWidth',1.5)
%         title('Chord')
%         legend('Storm','Hannah','Location','northeast')
%         xlabel('r (m)')
%         ylabel('Chord (m)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (c/63 - Hchord)./((c/63 + Hchord)/2),'LineWidth',1.5)
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
%     plot(r,sqrt(Wm),'b','LineWidth',1.5)
%         hold on
%     plot(Hr*63,HW,'r','LineWidth',1.5)
%         title('U_{rel}')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('r (m)')
%         ylabel('U_{rel} (m/s)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (sqrt(Wm) - HW)./((sqrt(Wm) + HW)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])

% figure;
% subplot(2,2,[1 2])
%     plot(r,CL,'b','LineWidth',1.5)
%         hold on
%     plot(Hr*63,test(1,:),'r','LineWidth',1.5)
%         title('Lift Coefficient')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('r (m)')
%         ylabel('CL (-)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (CL - test(1,:))./((CL + test(1,:))/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
% 
% figure;
% subplot(2,2,[1 2])
%     plot(r,CD,'b','LineWidth',1.5)
%         hold on
%     plot(Hr*63,Hcd,'r','LineWidth',1.5)
%         title('Drag Coefficient')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('r (m)')
%         ylabel('CD (-)')
%         xlim([min(r) max(r)])
% subplot(2,2,[3 4])
%     plot(r, 100 * (CD - Hcd)./((CD + Hcd)/2),'LineWidth',1.5)
%         xlabel('r (m)')
%         ylabel('Percent Difference (%)')
%         xlim([min(r) max(r)])
% 
% figure;
% subplot(2,2,[1 2])
%     plot(speed,sweep/1e6,'b','LineWidth',1.5)
%         hold on
%     plot(speed,Hannah,'r','LineWidth',1.5)
%         title('Power')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('Speed Shear Exponent (\alpha)')
%         ylabel('Power (MW)')
%         xlim([min(speed) max(speed)])
% subplot(2,2,[3 4])
%     plot(speed, 100 * (sweep/1e6 - Hannah)./((sweep/1e6 + Hannah)/2),'LineWidth',1.5)
%         xlabel('Speed Shear Exponent (\alpha)')
%         ylabel('Percent Difference (%)')
%         xlim([min(speed) max(speed)])

% figure;
% subplot(2,2,[1 2])
%     plot(1:length(r),r/63,'b','LineWidth',1.5)
%         hold on
%     plot(1:length(Hr),Hr,'r','LineWidth',1.5)
%         title('Radial Positions')
%         legend('Storm','Hannah','Location','northwest')
%         xlabel('Radial Position Index (-)')
%         ylabel('Radial Position (m)')
%         xlim([1 63])
% subplot(2,2,[3 4])
%     plot(1:length(Hr), 100 * (r/63 - Hr)./((r/63 + Hr)/2),'LineWidth',1.5)
%         xlabel('Radial Position Index (-)')
%         ylabel('Percent Difference (%)')
%         xlim([1 63])


