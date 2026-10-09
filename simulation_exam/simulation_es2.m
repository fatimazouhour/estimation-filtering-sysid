clear all
close all
clc

load data2.mat

A=[2.7258 , -1.2424 , 0.7577;
    2   ,0,  0;
    0  ,0.5 ,0];

B=[0.026,0,0]';
C=[1  -0.45  -0.05];
D=[0];

Bv1= 0.05*[0.026 0 0]';

V1= Bv1*Bv1';
V2=0.0004;
V12= 0;

P1=eye(3);
rng('default');

[n,nn]=size(A);

N=length(u);

x_avg=[0 0 0]';
x0= [50;100;50];

x(:,1) = x0;
for t=1:N
    v1(:,t)= mvnrnd(zeros(1,n),Bv1*Bv1');
    v2(t)= sqrt(V2)*randn;
    
    x(:,t+1)= A*x(:,t)+B*u(t)+v1(:,t);
    y(t) = C*x(:,t)+v2(t);
    
end

covv1= cov(v1'),V1;
covv2= cov(v2'),V2;

% ok
xhss(:,1) = x_avg;

S= ss(A,[B,eye(n)],C,[D,zeros(1,n)],1);

[kalman_pred,kbar,Pbar,kbar0]= kalman(S,V1,V2,0);

for t= 1:N 
    yhss(t)= C*xhss(:,t);
    ess(t) = y(t)-yhss(t);
    
    xhss(:,t+1)= A*xhss(:,t)+B*u(t)+kbar*ess(t);
    
    xfss(:,t)= xhss(:,t)+kbar0*ess(t);
    yfss(t)= C*xfss(:,t);
end


%---------kalman predictor----------------------------
xh(:,1)= x_avg;
P{1}= eye(3);

for t=1:N
   K0{t}=P{t}*C'*inv(C*P{t}*C'+V2);
   P0{t}= (eye(n)-K0{t}*C)*P{t};
   P{t+1}= A*P0{t}*A'+V1; 
   
   yh(t)= C*xh(:,t);
   e(t)= y(t)-C*xh(:,t);
   
   xf(:,t)= xh(:,t)+K0{t}*e(t);
   xh(:,t+1)= A*xf(:,t)+B*u(t);
  
end

K0= K0{t}
kbar

%-------RMSE----------------------------
N0_vector=[0  50 100];

for idx= 1:length(N0_vector)
   N0 =  N0_vector(idx);
   for k=1:n
   RMSE_xfss(k,idx)= norm(x(k,N0+1:N)-xfss(k,N0+1:N))/sqrt(N-N0);
   RMSE_xh(k,idx)= norm(x(k,N0+1:N)-xh(k,N0+1:N))/sqrt(N-N0);
   end
   RMSE_yfss(idx)= norm(y(N0+1:N)-yfss(N0+1:N))/sqrt(N-N0);
  % RMSE_yh(idx)= norm(y(N0+1:N)-yh(N0+1:N))/sqrt(N-N0);
end

RMSE_xfss, RMSE_xh,RMSE_yfss,%RMSE_yh