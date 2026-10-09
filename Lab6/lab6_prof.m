% example given by the prof
clear all
close all
clc

load heater.mat

figure
plot(ue)
figure
plot(ye)

%we remove the mean values

u1= ue-mean(ue);
y1= ye-mean(ye);


 u2= uv-mean(uv);
 y2= yv-mean(yv);
% we validate visually
figure, plot(u1), figure, plot(y1)
figure, plot(u2) , figure, plot(y2)
 
%% NN def
% 6 neurons --> 6*H  usually greater than 5 
% unique layer  output --> L;
NetDef= ['HHHHHH';'L-----'];

%% identification of NNARX ( na,nb,nk)
% perform 5 attempts and select the best 
na=1;
nb=na;
nk=1;

trparms=settrain;
trparms= settrain(trparms,'maxiter',500);
[W1,W2] = nnarx(NetDef,[na,nb,nk],[],[],trparms,y1',u1');
figure, drawnet(W1,W2) 

%% now we do the validation of the NNARX

alpha = W2(1:end-1);
alpha_0= W2(end);
beta=W1(:,1:end) ;
beta_0=W1(:,end) ;


tmin= max([na+1,nk+nb]);
tmax=length(y2);
for t=tmin:tmax
    phi= [y2(t-1:-1:t-na);u2(t-nk:-1:t-nk-nb+1)];
    yp(t) = alpha*tanh(beta*phi+beta_0)+alpha_0;
end


figure, plot(1:tmax,y2,'g',1:tmax,yp,'r-.')
N0=10;
na,nb,nk
Nv=tmax;
RMSE = sqrt(1/(Nv-N0))*norm(y2(N0+1:end)-yp(N0+1:end))
n = na+nb;
MDL= n*log(Nv-N0) / ( Nv-N0)+log(RMSE^2)



