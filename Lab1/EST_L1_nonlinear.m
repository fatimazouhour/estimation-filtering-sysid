clear all
close all
clc

 % same ex but estimating it all  includeing the non linear part 
load sensor.mat

format short e % law ma 3mlt hek ma kent byn el nbs bl theta_hat

N=length(z);
phi = [Vz.^3,Vz.^2,Vz,Vz.^0];

theta_hat = phi\z;
Vz0= linspace(-10,10,1000);

z_pol = polyval(theta_hat,Vz0)

figure, plot(z,Vz,'o',z_pol,Vz0,'--'),
axis([min(z),max(z),-10,10])

grid on