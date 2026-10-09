clear all
close all
clc

load data.mat

A=[0.96 0.5 0.27 0.28;
    -0.125 0.96 -0.08 -0.07;
    0 0 0.85 0.97;
    0 0 0 0.99];

B=[1 -1 2 1]';

C= [0 2 0 0];
D=[0];

rng('default');

Bv1= sqrt(15)*[0.5, 0,0,1]';
V2= 2000;
V1= Bv1*Bv1';
V12= 0;

[n,nn]= size(A);
N= length(u);

P{1}= 0.5*eye(n);

x(:,1)= [30 40 -70 -10]';

for t=1:N
   v1(:,t)= mvnrnd(zeros(1,n),Bv1*Bv1');
   v2(t)= sqrt(V2)*randn;
   
   x(:,t+1)= A*x(:,t)+B*u(t)+v1(:,t);
   y(t)= C*x(:,t)+v2(t);
   
end

% dynamic filter 

xh(:,1)= zeros(n,1);

for t=1:N
    yh(t)= C*xh(:,t);
    e(t)= y(t)-yh(t);
    
    K0{t}= P{t}*C'*inv(C*P{t}*C'+V2);
    K{t}= (A*P{t}*C')*inv(C*P{t}*C'+V2);
    P{t+1}= A*P{t}*A'+V1-K{t}*(C*P{t}*C'+V2)*K{t}';
    
    xh(:,t+1)=A*xh(:,t)+B*u(t)+K{t}*e(t);
    xf(:,t)= xh(:,t)+K0{t}*e(t);
    
    yf(t)= C*xf(:,t);

end
K0{N},K{N}

% steady state

xhss(:,1)= zeros(n,1);
S=ss(A,[B,eye(n)],C,[D,zeros(1,n)],1);
[kalman_pred,Kbar,Pbar,K0bar]=kalman(S,V1,V2,0);


for t=1:N
    yhss(t)= C*xhss(:,t);

    e(t)= y(t)-yhss(t);
    
    xhss(:,t+1)= A*xhss(:,t)+B*u(t)+Kbar*e(t);
    xfss(:,t)= xhss(:,t)+K0bar*e(t);
    
    yfss(t)= C*xfss(:,t);
    
end

K0bar,Kbar

% dynamix in predictor correctorr form

xhc(:,1)=zeros(n,1); 
Pc{1}= 0.5*eye(4);

for t=1:N
   K0c{t}= Pc{t}*C'*inv(C*Pc{t}*C'+V2);
   P0{t}= (eye(n)-K0c{t}*C)*Pc{t};
   Pc{t+1}= A*P0{t}*A'+V1;
   
   e(t)= y(t)-C*xhc(:,t);
   
   xfc(:,t)= xhc(:,t)+K0c{t}*e(t);
   xhc(:,t+1)= A*xfc(:,t)+B*u(t);
    
end

norm(xh-xhc,inf)
norm(xf-xfc,inf)
% to check the quality we see rmse 

N0_vector= [0,20,100];

for idx= 1:length(N0_vector)
    N0= N0_vector(idx);
    for k= 1:n
        RMSE_xh(k,idx)= norm(x(k,N0+1:N)-xh(k,N0+1:N))*sqrt(1/(N-N0));
         RMSE_xf(k,idx)= norm(x(k,N0+1:N)-xf(k,N0+1:N))*sqrt(1/(N-N0));
    end
    RMSE_yh(idx)= norm(y(N0+1:N)-yh(N0+1:N))*sqrt(1/(N-N0));
     RMSE_yf(idx)= norm(y(N0+1:N)-yf(N0+1:N))*sqrt(1/(N-N0));
end

RMSE_yh, RMSE_xh,RMSE_xf,RMSE_yf

for idx= 1:length(N0_vector)
    N0= N0_vector(idx);
    for k= 1:n
        RMSE_xhss(k,idx)= norm(x(k,N0+1:N)-xhss(k,N0+1:N))*sqrt(1/(N-N0));
         RMSE_xfss(k,idx)= norm(x(k,N0+1:N)-xfss(k,N0+1:N))*sqrt(1/(N-N0));
    end
    RMSE_yhss(idx)= norm(y(N0+1:N)-yhss(N0+1:N))*sqrt(1/(N-N0));
     RMSE_yfss(idx)= norm(y(N0+1:N)-yfss(N0+1:N))*sqrt(1/(N-N0));
end

RMSE_xhss,RMSE_yhss,RMSE_xfss,RMSE_yfss

