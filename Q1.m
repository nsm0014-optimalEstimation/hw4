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
tf = 25;
dt = 0.005;
time = ts:dt:tf;

theta_ddot = zeros(1,length(time));
theta_dot = zeros(1,length(time));
theta_dot_3 = zeros(1,length(time));
theta = zeros(1,length(time));
f_dist = 5 + gaussianDistFCN([1,length(time)],sqrt(2),0);
y = zeros(1,length(time));

for i = 2:length(time)

    theta_ddot(i) = tau/Jp + f_dist(i)*L*theta(i-1)/Jp - m*g*L*theta(i-1)/Jp - b*theta_dot(i-1)^3/Jp;
    theta_dot(i) = theta_dot(i-1) + theta_ddot(i)*dt;
    theta(i) = theta(i-1) + theta_dot(i)*dt;
    theta_dot_3(i) = theta_dot(i)^3;

    y(i) = theta(i) + gaussianDistFCN([1 1],sqrt(deg2rad(1)),0);

end

plot(time,y)
figure
plot(time,theta_dot_3)