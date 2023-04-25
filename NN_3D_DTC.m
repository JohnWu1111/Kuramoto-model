% a demo of nearest neighbor 3D Kuramoto model

clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

%% parameter
L = 10;
K = 0.5;
delta = 0.5;
w = 0.1;
omega = randn(L,L,L);
theta0 = 2*pi*rand(L,L,L);
T = 100;
dt = 1e-2;
t = 0:dt:T;
nt = length(t);
order = zeros(1,nt);

%% time evolution
theta = theta0;
order(1) = abs(sum(exp(1i*theta),"all"))/L^3;
for i = 2:nt
    theta = myrunge(theta,dt,omega,(K+delta*cos(2*pi*w*t(i))));
    order(i) = abs(sum(exp(1i*theta),"all"))/L^3;
end

%% analysis and plot
order_mean = mean(order(:,floor(nt/2):end),2);

figure;
plot(t,order)


toc;

%% functions
function y = myrunge(x,dt,omega,fact)
c1 = coeff(x,omega,fact);
c2 = coeff(x+c1*dt/2,omega,fact);
c3 = coeff(x+c2*dt/2,omega,fact);
c4 = coeff(x+c3*dt,omega,fact);
y = x + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = coeff(x,omega,fact)
y = omega - fact*(sin(x-circshift(x,1)) + sin(x-circshift(x,-1)) ...
    + sin(x-circshift(x,1,2)) + sin(x-circshift(x,-1,2)) ...
    + sin(x-circshift(x,1,3)) + sin(x-circshift(x,-1,3)));
end