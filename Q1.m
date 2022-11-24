%% Formatting
clc
clear
close all
format shortg

%% Begin Question 1
Jp = 2.5;
m  = 1.6;
L = 1;
b = 1.25;
tau = 12;
g = 9.81;

ts = 0;
tf = 50;
dt = 0.01;
time = ts:dt:tf;

theta_ddot = zeros(1,length(time));
theta_dot = zeros(1,length(time));
theta_dot_3 = zeros(1,length(time));
theta = zeros(1,length(time));
f_dist = 5 + gaussianDistFCN([1,length(time)],sqrt(2),0);
y = zeros(1,length(time));

A_discrete = c2dm([-b/Jp, -(m*g*L/Jp) + L*f_dist(1)/Jp; 1,0],[1/Jp;0],[0 1], [0], 0.001)
A_expm = expm(A_discrete)
for i = 2:length(time)

    theta_ddot(i) = tau/Jp + f_dist(i)*L*theta(i-1)/Jp - m*g*L*theta(i-1)/Jp - b*theta_dot(i-1)^3/Jp;
    theta_dot(i) = theta_dot(i-1) + theta_ddot(i)*dt;
    theta(i) = theta(i-1) + theta_dot(i)*dt;
    theta_dot_3(i) = theta_dot(i)^3;

    y(i) = theta(i) + gaussianDistFCN([1 1],sqrt(deg2rad(1)),0);

end


modelFig = figure('Position',[250 100 1000 600]);
hold on
plot(time,y)
plot(time,theta,'LineWidth',2)
xlabel('Time [s]')
ylabel('Deflection [deg]')
title('Modeled Pendelum with a Disturbance Force')
legend('Measurement','Model')
ylim([0 2])
ax = gca;
ax.FontSize = 18;
saveas(modelFig,'Q1model.png')

syms theta_dot theta f_dist
theta_ddot = tau/Jp + f_dist*L*theta/Jp - m*g*L*theta/Jp - b*theta_dot^3/Jp
int(theta_ddot)
pretty(ans)
pretty(theta_ddot)
