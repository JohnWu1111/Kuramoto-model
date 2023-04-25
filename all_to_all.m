% a demo of all-to-all Kuramoto model

clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed)

%% parameter
L = 1000;
K_all = 0.3:0.3:3;
nK = length(K_all);
omega = randn(L,1);
theta0 = 2*pi*rand(L,1);
T = 100;
dt = 1e-2;
t = 0:dt:T;
nt = length(t);
order = zeros(nK,nt);

%% time evolition
for n = 1:nK
    K = K_all(n);
    theta = theta0;
    order(n,1) = sum(exp(1i*theta))/L;
    for i = 2:nt
        theta = myrunge(theta,K*order(n,i-1),dt,omega);
        order(n,i) = sum(exp(1i*theta))/L;
    end    
end

%% analysis and plot
order_mean = mean(abs(order(:,floor(nt/2):end)),2);

figure;
subplot(2,1,1)
plot(t,abs(order))
subplot(2,1,2)
plot(K_all,order_mean)

toc;

%% functions
function y = myrunge(x,r,dt,omega)
c1 = coeff(x,r,omega);
c2 = coeff(x+c1*dt/2,r,omega);
c3 = coeff(x+c2*dt/2,r,omega);
c4 = coeff(x+c3*dt,r,omega);
y = x + dt*(c1+2*c2+2*c3+c4)/6;
end

function y = coeff(x,r,omega)
y = omega + abs(r)*sin(angle(r)-x);
end