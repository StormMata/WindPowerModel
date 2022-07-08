%% TSR Sensitivity

sp = ('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Research/WindPowerModel/5 MW Turbine/DShear Sensitivity');

addpath(sp);

load('Turbine.mat')

%% change I

% alp = -0.1:0.1:0.7;
% 
% % for i = 1:length(alp)
% %     I.SpeedShear(:,i) = 7.75 * (I.Heights./T.Hub).^alp(i);
% % end
% 
% I.SpeedShear = 9 * ones(length(I.Heights),size(I.SpeedShear,2));
% I.DirShear   = zeros(length(I.Heights),size(I.SpeedShear,2));
% I.Omega      = 9*7.5/63 * ones(1,size(I.SpeedShear,2));

%% Generate Test Range - Linear Direction Shear

    DType = ('Linear');

    U = 9;
    
    TSR = 10 - 5 * (0.6847/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/DSsweep_%s.mat',sp,DType),'I')

%% Generate Test Range - Tangent Direction Shear

    DType = ('Tangent');

    U = 9;
    
    TSR = 10 - 5 * (0.6847/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    x = linspace(-pi/2.5,pi/2.5,126);
    y = tan(x);
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
    
            ytemp                 = (y/max(y)*T.R*direc(i,j))';
            ytemp                 = [ytemp(end)*ones(47,1); flip(ytemp); ytemp(1)*ones(27,1)];
            I.DirShear(:,index)   = ytemp - ytemp(111);
    
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/DSsweep_%s.mat',sp,DType),'I')

%% Generate Test Range - Arc Direction Shear

    DType = ('Arc');

    U = 9;
    
    TSR = 10 - 5 * (0.6847/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';

    r = 100;

    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);

            if direc(i,j) >= 0.6
                r = 100;
            elseif direc(i,j) == 0.4 || direc(i,j) == 0.5
                r = 110;
            elseif direc(i,j) == 0.3
                r = 105;
            elseif direc(i,j) == 0.2
                r = 120;
            elseif direc(i,j) == 0.0 || direc(i,j) == 0.1
                r = 110;
            elseif direc(i,j) == -0.1
                r = 70;
            end
    
            p = [-T.R*direc(i,j) T.R*direc(i,j); 0 126];
    
            a = sym('a',[2,1],'real');
            eqs = [1,1]*(p - repmat(a(:),1,2)).^2 - r^2;
            sol = vpasolve(eqs,a);
            ss = struct2cell(sol);
            xy = double([ss{:}]);
            
            v = xy(1,:);
            p1 = p - v(:);
            alp = atand(p1(2,:)./p1(1,:));
            alp = alp + 180*(p1(1,:) < 0 & p1(2,:) > 0) - 180*(p1(1,:) < 0 & p1(2,:) < 0);
            asort = sort(alp);
            phi = linspace(asort(1),asort(2),126)';

            y = r*cosd(phi) + v(1);
    
            I.DirShear(:,index)   = [y(end)*ones(47,1); flip(y); y(1)*ones(27,1)];
    
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;

        end
    end

    save(sprintf('%s/DSsweep_%s.mat',sp,DType),'I')

% % % ------
% p = [-T.R*-0.3 T.R*0.3; 0 126];
% r = 105;
% 
% a = sym('a',[2,1],'real');
% eqs = [1,1]*(p - repmat(a(:),1,2)).^2 - r^2;
% sol = vpasolve(eqs,a);
% ss = struct2cell(sol);
% xy = double([ss{:}]);
% 
% v = xy(1,:);
% p1 = p - v(:);
% alp = atand(p1(2,:)./p1(1,:));
% alp = alp + 180*(p1(1,:) < 0 & p1(2,:) > 0) - 180*(p1(1,:) < 0 & p1(2,:) < 0);
% asort = sort(alp);
% phi = linspace(asort(1),asort(2),126)';
% hold on
% plot(r*cosd(phi) + v(1),r*sind(phi) + v(2),'-b',v(1),v(2),'ok')
% grid on
% 
% y = r*cosd(phi) + v(1);
% 
% y = [y(end)*ones(47,1); flip(y); y(1)*ones(27,1)] .* ones(200,10);
% 
% I.DirShear(:,41:50) = y;

%% Generate Test Range - Arc Direction Shear

x = -10:0.1:10;
f = polyfit([-5 0 5],[0 T.R 126],6)
% y = f(1)*x.^2 + f(2)*x.^1 + f(3)*x.^0;
% y = f(1)*x.^3 + f(2)*x.^2 + f(3)*x.^1 + f(4)*x.^0;
% y = f(1)*x.^4 + f(2)*x.^3 + f(3)*x.^2 + f(4)*x.^1 + f(5)*x.^0;
% y = f(1)*x.^5 + f(2)*x.^4 + f(3)*x.^3 + f(4)*x.^2 + f(5)*x.^1 + f(6)*x.^0;
y = f(1)*x.^6 + f(2)*x.^5 + f(3)*x.^4 + f(4)*x.^3 + f(5)*x.^2 + f(6)*x.^1 + f(7)*x.^0;
plot(x,y)
hold on
xline(-5)
xline(5)
yline(0)
yline(126)
scatter(0,T.R,20,"red","filled")
xlim([-10 10])
ylim([-200 300])

%% Generate Test Range - Random Direction Shear

[speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
index          = 1;

r = 10*(rand(200,8)-0.5);

for i = 1:size(direc,1)
    for j = 1:size(speed,2)
        I.SpeedShear(:,index) = 9 * (I.Heights./T.Hub).^speed(i,j);

        ytemp                 = direc(i,j)*I.Heights + r(:,i)*direc(i,j);

        I.DirShear(:,index)   = ytemp - ytemp(111);

        I.Omega(index)        = 9*7.5/63;
        index = index + 1;
    end
end

plot(I.DirShear(:,1),I.Heights,'LineWidth',2)
hold on
yline(T.Hub+T.R)
yline(T.Hub-T.R)
xline(nonzeros(I.DirShear(:,1).*(I.Heights==(T.Hub+T.R))))
xline(nonzeros(I.DirShear(:,1).*(I.Heights==(T.Hub-T.R))))
xlabel('Wind Direction (deg)')
ylabel('z (m)')

%% Generate Test Range - Random Direction Shear Nonlinear

[speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
index          = 1;

r = 10*(rand(200,8)-0.5);

for i = 1:size(direc,1)
    for j = 1:size(speed,2)
        I.SpeedShear(:,index) = 9 * (I.Heights./T.Hub).^speed(i,j);

        ytemp                 = direc(i,j)*I.Heights + r(:,i)*direc(i,j);

        I.DirShear(:,index)   = ytemp - ytemp(111);

        I.Omega(index)        = 9*7.5/63;
        index = index + 1;
    end
end

plot(I.DirShear(:,1),I.Heights,'LineWidth',2)
hold on
yline(T.Hub+T.R)
yline(T.Hub-T.R)
xline(nonzeros(I.DirShear(:,1).*(I.Heights==(T.Hub+T.R))))
xline(nonzeros(I.DirShear(:,1).*(I.Heights==(T.Hub-T.R))))
xlabel('Wind Direction (deg)')
ylabel('z (m)')