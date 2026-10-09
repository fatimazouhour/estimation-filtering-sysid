clear all
close all
clc

format compact, format short

% same as before but using set membership theory 

load sensor.mat

i1 =14;
i2= 36;

z_linear = z(i1:i2);
Vz_linear = Vz(i1:i2);

N= length(z_linear);

phi = [Vz_linear, ones(N,1)];

theta_hat = phi\z_linear;

kt_hat= 1/ theta_hat(1)
V0_hat = -theta_hat(2)/theta_hat(1)

Vz0= linspace(-10,10,1000);

z_hat = theta_hat(1)*Vz0+theta_hat(2);


% i can plot these results
figure, plot(z,Vz,'o', z_hat,Vz0,'-', ...
             z_linear(1)*[1,1],[-10,10],'-', z_linear(end)*[1,1],[-10,10],'-'),

grid on
 
% now i want to evaluate the uncertainty of the values 
% EUI_inf

A= pinv(phi);
e= 0.5*1e-3;

format shorte
for i= 1:length(theta_hat)
    theta_hat_m(i,1)= A(i,:)*(z_linear-e*sign(A(i,:))');
end

theta_hat_M = 2*theta_hat-theta_hat_m;

EUI_inf= [theta_hat_m,theta_hat_M]
format short
% fina nst3ml loop bs 3aj2a no need really
% now for the kt and V0 

kt_interval=[1/ theta_hat_M(1), 1/theta_hat_m(1)]
V0_interval = [-theta_hat_M(2)/theta_hat_m(1),-theta_hat_m(2)/theta_hat_M(1)]

z_min = min([theta_hat_M(1)*Vz0;theta_hat_m(1)*Vz0])+theta_hat_m(2);
z_max = max([theta_hat_m(1)*Vz0;theta_hat_M(1)*Vz0])+theta_hat_M(2);

% i can plot these
figure
plot(z,Vz,'o',z_hat,Vz0,'-',z_min,Vz0,'--',z_max,Vz0,'--',z_linear(1)*[1,1],[-10,10],'-', z_linear(end)*[1,1],[-10,10],'-')
grid on

% i can also use the PUI to estimate 

M=[phi;-phi];
b = [z_linear+e;-z_linear+e];
format shorte

for i = 1:2
    c = zeros(2,1); c(i) = 1;
    PUI(i,1) = c' * linprog( c, M, b);     % min
    PUI(i,2) = c' * linprog(-c, M, b);     % max
end

% now to get kt and Vo
format short
kt_pui = [1/PUI(1,2),1/PUI(1,1)];
Vo_pui = [-PUI(2,2)/PUI(1,1),-PUI(2,1)/PUI(1,2)];

% we can get the central estimates
format short 
pui_central = [(PUI(1,1)+PUI(1,2))/2;(PUI(2,1)+PUI(2,2))/2]
kt_central = 1/pui_central(1)
Vo_central= -pui_central(2)/pui_central(1)

% hence we can get z_max and z_min central 
z_central = Vz0/kt_central -Vo_central/kt_central;
z_min_central = min([PUI(1,2)*Vz0;PUI(1,1)*Vz0])+PUI(2,1);
z_max_central = max([PUI(1,2)*Vz0;PUI(1,1)*Vz0])+PUI(2,2);

figure,
plot(z,Vz,'o',z_central,Vz0,'-',z_min_central,Vz0,'--',z_max_central,Vz0,'--',z_linear(1)*[1,1],[-10,10], ...
    '-',z_linear(end)*[1,1],[-10,10],'-')
grid on

% i can get as well the feasable parameter set FPS

% let n=16 be the number of the polytype 
% IMP  ASK PROF ABOUT IT 

n = 16;

for k = 1:n
   c= [sin(2*pi*k/n),cos(2*pi*k/n)];
   theta_vect(:,k)= linprog(c,M,b);
end

k=convhull(theta_vect(1,:),theta_vect(2,:));

% then plot
