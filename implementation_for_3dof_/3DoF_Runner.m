function DoF_Runner()
% Mahmoud Muhammad Yahaya's 3DOF  model implementation solve via NLS algo.
% mahmoudpd@gmail.com
clear;
t0 = 0;
tf = 10; % maximum sec.
g = 0.05; % sampling period
Kmax = (tf - t0)/g; % 

theta0 = [0, pi/3, pi/2]; % the angle is measeured in radians
thetamah = theta0;

tic;
for k = 1:Kmax
    tk = k*g;
    rd = [1.5 + 0.2*sin(1*tk), sqrt(3)/2 + 0.2*sin(2*tk + pi/7)]'; % infinity
    thetak = SDMSCfor3DOF(rd, thetamah(end,:)'); % calling the sdmsc algorithm 
    thetamah = [thetamah; thetak']; % a matrix of joint angular vectors
end
toc
%disp([num2str(toc)]);

% data = load('pseudo_data.m');

x0 = zeros(Kmax+1,3);
x1 = [cos(thetamah(:,1)), sin(thetamah(:,1))];
x2 = [cos(thetamah(:,1)) + cos(thetamah(:,1) + thetamah(:,2)), sin(thetamah(:,1)) + sin(thetamah(:,1) + thetamah(:,2))];
x3 = [cos(thetamah(:,1)) + cos(thetamah(:,1) + thetamah(:,2)) + cos(thetamah(:,1) + thetamah(:,2) + thetamah(:,3)), sin(thetamah(:,1)) + sin(thetamah(:,1) + thetamah(:,2)) + sin(thetamah(:,1) + thetamah(:,2)+ thetamah(:,3))];


% simulate the robot arm  
for k = 1:Kmax
    plot([x0(k,1); x1(k,1)],[x0(k,2); x1(k,2)],'r');
    grid on
    hold on
    plot([x1(k,1); x2(k,1)],[x1(k,2); x2(k,2)],'g');
    grid on
    hold on
    plot([x2(k,1); x3(k,1)],[x2(k,2); x3(k,2)],'b');
    grid on
end

% plot that draw the Lassajous curve at the end-effector
rdz = [];
for k = 0:Kmax
    tk = k*g;
    rd = [1.5 + 0.2*sin(1*tk), sqrt(3)/2 + 0.2*sin(2*tk + pi/7)]; % infinity
    rdz = [rdz; rd];
end
plot(rdz(:,1),rdz(:,2))
grid on

% show-case the path tracking of the Lissajous by the solve model
figure
plot(rdz(:,1),rdz(:,2),'ob', x3(:,1),x3(:,2),'-.r');
legend('Desired path','Actual trajectory')
grid on

figure
semilogy(abs(x3(:,1) - rdz(:,1)), 'MarkerSize',3,'LineWidth',2)
grid on
axis([0,200,0,1])
legend('\rm SDMSC: along x-coord')
grid on

figure
semilogy(abs(x3(:,2) - rdz(:,2)), 'MarkerSize',3,'LineWidth',2)
grid on
axis([0,200,0,1])
legend('\rm SDMSC: along y-coord ')
grid on



