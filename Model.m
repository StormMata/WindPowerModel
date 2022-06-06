clearvars -except AeroDyn A
yaw     = 0;
azimuth = linspace(0,2*pi,360);
U       = 11.4*(1-1/3);
r       = AeroDyn.r+1.5;
Omega   = (12.1*2*pi)/(60);
rho     = 1.225;

for j = 1:length(azimuth)

    for i = 1:length(r)
    
        W(i)   = (U * cos(yaw * sin(azimuth(j))) * cos(yaw * cos(azimuth(j))))^2 + ...
                 (Omega * r(i) - U * cos(yaw * sin(azimuth(j))) * sin(yaw * cos(azimuth(j))))^2;
    
        Phi(i) = atan((U * cos(yaw * sin(azimuth(j))) * cos(yaw * cos(azimuth(j))))/ ...
                      (Omega * r(i) - U * cos(yaw * sin(azimuth(j))) * sin(yaw * cos(azimuth(j)))));
    
        c      = AeroDyn.Chord(i);
    
        alpha  = rad2deg(Phi(i)) - AeroDyn.Twist(i);
    
%         row    = find(A.(sprintf('%s',AeroDyn.AeroID(i))).Alpha == alpha);
    
%         if isempty(row)
%     
            p1 = find(A.(sprintf('%s',AeroDyn.AeroID(i))).Alpha > alpha, 1, 'first');
    
            CL = sum(A.(sprintf('%s',AeroDyn.AeroID(i))).CL(p1-1:p1))/2;
            CD = sum(A.(sprintf('%s',AeroDyn.AeroID(i))).CD(p1-1:p1))/2;
%     
%         else
    
%             CL = A.(sprintf('%s',AeroDyn.AeroID(i))).CL(row);
%             CD = A.(sprintf('%s',AeroDyn.AeroID(i))).CD(row);
    
%         end
    
        dQ(i)  = 3 * 0.5 * rho * c * W(i) * r(i) * (CL * sin(Phi(i)) - CD * cos(Phi(i)));
    
    end

dQ = trapz(r,dQ);

dP(j) = Omega * dQ;

end

trapz(azimuth,dP)/(2*pi)

% A.(sprintf('%s',AeroDyn.AeroID(i))).Alpha
% A.NACA64A17.Alpha

% table2array(AeroDyn((AeroDyn.r==r(i)),3));

% y = interp1(A.DU35A17.Alpha(),A.DU35A17.CL(),-180:180,"spline");