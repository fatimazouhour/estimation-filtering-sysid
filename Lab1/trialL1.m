clear all
close all
clc

load sensor.mat

idx1=11;
idx2=36;

figure,
plot(z,Vz,'o')
hold on, grid on
xlabel('z'),ylabel('Vz')


z_linear=z(idx1:idx2);
Vz_linear= Vz(idx1:idx2);

N=length(Vz_linear);

phi = [Vz_linear ones(N,1)];

theta= phi\z_linear;

kt_hat= 1/theta(1);
Vo_hat= theta(2)*kt_hat;

z_hat = Vo_hat*1/kt_hat+Vz_linear*1/kt_hat;

figure, plot(z_hat,Vz_linear,'o'), xlabel('z_linear'),ylabel('Vz_linear')
