clear all
close all
clc
format compact

load sensor.mat  % i have the parameters Vz & z

figure(1)
plot(z,Vz,'Marker','o')
hold on, grid on
xlabel('z')
ylabel('Vz')

% notice that we have linearity b/n 0.01 & 0.036
% so i want to work only on that vector for now
% defining: 
idx1= 11;
idx2= 36;

z_linear= z(idx1:idx2);
Vz_linear = Vz(idx1:idx2);

N = length(z_linear);

phi = [Vz_linear ones(N,1)];

% LEAST SQUARE ESTIMATE
theta_hat= phi\z_linear;
kt_hat= 1/theta_hat(1)
Vo_hat = -theta_hat(2)/theta_hat(1)

% now plot the above result 
Vz0 = linspace(-10,10,10000);
z_hat = theta_hat(1)*Vz0+theta_hat(2); % nfs el eq el aseseye bs estimation

figure(2)
plot(z,Vz,'o',z_hat,Vz0,'-',z_linear(1)*[1,1],[-10,10],'-', ...
    z_linear(end)*[1,1],[-10,10],'-')
hold on
axis([min(z),max(z),-10,10])
hold on,grid on

% assuming P=95%
k= 2;
k_theta= 2;
noise_max= 0.5*1e-3;

sigma_e= noise_max/k;

sigma_e= sigma_e^2*eye(N) % sum matrix sigma
%convariance matrix 

covar_theta= inv(phi'*inv(sigma_e)*phi);

sigma_theta= sqrt(diag(covar_theta));
% IMP notice we have to multiply by 2 the sigma

theta_max= theta_hat+sigma_theta*2;
theta_min= theta_hat-sigma_theta*2;

% now for the kt

kt_min= 1/theta_max(1)
kt_max= 1/theta_min(1)

% Vo
Vo_max= -theta_max(2)/theta_min(1)
Vo_min= -theta_min(2)/ theta_max(1)

% we can plot the estimation
z_min=min([theta_max(1)*Vz0;theta_min(1)*Vz0])+theta_min(2);
z_max=max([theta_max(1)*Vz0;theta_min(1)*Vz0])+theta_max(2);

figure(3)
plot(z,Vz,'o',z_hat,Vz0,'-',z_min,Vz0,'--',z_max,Vz0,'-.', ...
    z_linear(1)*[1,1],[-10,10],'-',z_linear(end)*[1,1],[-10,10],'-')
hold on, grid on


%% Trial lab1 
clear all
close all
clc

format compact

load sensor.mat

i_11=11;
i_36 = 36;

Vz = Vz(i_11:i_36);
z = z(i_11:i_36); % since we want to only approximate the linear line for now

phi = [Vz, ones(26,1)]; 

% in our case y =z 
% where y= phi*theta +e

theta_hat= phi\z;

kt_hat = 1/theta_hat(1);
Vo_hat = -theta_hat(2)*kt_hat;

% i have the estimations for kt and Vo

% but how much do i trust those estimates? 
% lets use the prior info that e=0.5mm max and work from that 

sigma_e= 0.25*1e-3;
covar_theta = sigma_e^2 * inv(phi'*phi);

sigma_theta_1= sqrt(covar_theta(1,1));
sigma_theta_2= sqrt(covar_theta(2,2));

% now we can build confidence parameters

delta_theta_1= 2*sigma_theta_1;
delta_theta_2= 2*sigma_theta_2;

theta_hat1_min= theta_hat(1)-delta_theta_1;
theta_hat1_max = theta_hat(1) +delta_theta_1;


theta_hat2_min= theta_hat(2)-delta_theta_2;
theta_hat2_max= theta_hat(2)+delta_theta_2;


% now since we got the conf intervals for theta_hat1 and theta_hat2
% we can get easily for Vo and kt

kt_hat_interval= [1/theta_hat1_max;1/theta_hat1_min];
Vo_hat_interval = [-theta_hat2_max/theta_hat1_min;-theta_hat2_min/theta_hat1_max];

% fine kmn a3ml calc lal z_hat
Vz0= linspace(-10,10,1000);
z_hat = theta_hat(1)*Vz0+theta_hat(2);
% we can plot the figure but first i should get zmin and zmax

z_min = theta_hat1_min*Vz0+theta_hat2_min
z_max= theta_hat1_max*Vz0 +theta_hat2_max


figure(3)
plot(z,Vz,'o',z_hat,Vz0,'-',z_min,Vz0,'--',z_max,Vz0,'-.', ...
    z(1)*[1,1],[-10,10],'-',z(end)*[1,1],[-10,10],'-') % check this one out 
grid on

% i can consider as well the data driven in stead of the priori knowledge
k= 2;
k_theta= 2;
noise_max= 0.5*1e-3;

N = length(z);
sigma_e= noise_max/k;
sigma_e= sigma_e^2*eye(N)
