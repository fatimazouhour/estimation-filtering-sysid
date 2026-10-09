
clear all
close all
clc

load data.mat
% matrices given 
A= [0.96 0.5 0.27 0.28;-0.125 0.96 -0.08 -0.07;
    0 0 0.85 0.97;0 0 0 0.99];
B= [1 -1 2 1]';
C= [0 2 0 0];
D=[0];

% white noises and their matrices 
B_v1= sqrt(15)*[0.5 0 0 1]';
V1=B_v1*B_v1';
V2=2000; % variance
V12=0; % they are not correlated 

%
N=length(u);
N0=[0,20,100];
P1= 0.5*eye(4);


rng('default'); % this makes it produce the same random numbers 

% LTI system initialization
x(:,1) = [30 40 -70 -10]';
[n,m]= size(A);
for t = 1:N
    v1(:,t) = mvnrnd(zeros(1,n),B_v1*B_v1')';
    v2(t) = sqrt(V2)*randn;
    x(:,t+1)= A*x(:,t)+B*u(t)+v1(:,t);
    y(t) = C*x(:,t) + D*u(t) +v2(t);
end

%% Plot LTI system states 
% T=1:N;
% for k= 1:n
%     figure, plot(T,x(k,T),'g'), title('state x')
% end
% 
% figure, plot(T,y(T),'g') , title('output y ')
% 
%% DYNAMIC KALMAN 1_STEP PREDICTOR
xh(:,1) = zeros(n,1);
P{1}=P1;

for t= 1:N
    yh(t)=C*xh(:,t);
    e(t)= y(t)-yh(t);
    K{t}=(A*P{t}*C'+V12)*inv(C*P{t}*C'+V2);
    xh(:,t+1)= A*xh(:,t) +B*u(t)+K{t}*e(t);
    P{t+1}= A*P{t}*A' + V1 - K{t}*(C*P{t}*C'+V2)*K{t}';    
end

%% DYNAMIC KALMAN FILTER 
for t=1:N
    K0{t}=P{t}*C'*inv(C*P{t}*C'+V2);
    % k{t} from before is the same
    xf(:,t)= xh(:,t)+K0{t}*e(t);
    yf(t)= C*xf(:,t);
end

K_N=K{N}
K0_N=K0{N}
%% now we need to get the Root mean square and compare:

for idx=1:length(N0)
    N0_idx= N0(idx);
    for k=1:n,
        RMSE_xh(k,idx) = norm(x(k,N0_idx+1:N)- xh(k,N0_idx+1:N))/sqrt(N-N0_idx);
        RMSE_xf(k,idx) = norm(x(k,N0_idx+1:N)- xf(k,N0_idx+1:N))/sqrt(N-N0_idx);
    end
    RMSE_yh(k,idx) = norm(y(N0_idx+1:N)- yh(N0_idx+1:N))/sqrt(N-N0_idx);
    RMSE_yf(k,idx) = norm(y(N0_idx+1:N)- yf(N0_idx+1:N))/sqrt(N-N0_idx);

end

N0(1),N0(2),N0(3)
RMSE_xh, RMSE_yh, RMSE_xf, RMSE_yf % t2kde menon!!

%% Plot it all together 
T=1:N;
for k=1:n
    figure, plot(T,x(k,1:N),'g',T,xh(k,1:N),'b-.',T,xf(k,1:N),'r--')
    legend('s','predictor k', 'filter f'),title(['State x_', num2str(k),'(t)'])
end

figure, plot(T,y(1:N),'g', T,yh(1:N),'b-.',T,yf(1:N),'r--')
legend('s','predictor k', 'filter f')

%% Now for the steady state kalman filter and predictor 
xh_inf(:,1)=zeros(n,1);
S= ss(A,[B, eye(n)], C,[D,zeros(1,n)],1);
[Kalman_p,Kbar,Pbar,K0bar]= kalman(S,V1,V2,0)

for t=1:N
    yh_inf(t)= C*xh_inf(:,t);
    e_inf(t)= y(t)- yh_inf(t);
    xh_inf(:,t+1)=A*xh_inf(:,t)+B*u(t)+Kbar*e_inf(t);
    
    % for the filter 
    xf_inf(:,t)= xh_inf(:,t)+ K0bar*e_inf(t);
    yf_inf(t)= C*xf_inf(:,t);
end

% now RMSE 
for idx= 1:length(N0)
    N0_idx= N0(idx);
    for k=1:n
        RMSE_xh_inf(k,idx)= norm(x(k,N0_idx+1:N)-xh_inf(k,N0_idx+1:N))/sqrt(N-N0_idx);
        RMSE_xf_inf(k,idx)= norm(x(k,N0_idx+1:N)-xf_inf(k,N0_idx+1:N))/sqrt(N-N0_idx);
    end
    RMSE_yh_inf(idx)= norm(y(N0_idx+1:N)-yh_inf(N0_idx+1:N))/sqrt(N-N0_idx);
    RMSE_yf_inf(idx)= norm(y(N0_idx+1:N)-yf_inf(N0_idx+1:N))/sqrt(N-N0_idx);
end
RMSE_xh_inf, RMSE_xf_inf,RMSE_yh_inf,RMSE_yf_inf


%% ALso do the plot here 
% ha a3melo skip b3den brj3lo

%% Dynamic predictor Kpc in predictor corrector form 
xh_pc(:,1) = zeros(n,1);
P_pc{1}= P1;

for t=1:N
    K0_pc{t}= P_pc{t}*C'*inv(C*P_pc{t}*C'+V2);
    P0_pc{t}= (eye(n)-K0_pc{t}*C)*P_pc{t};
    e_pc(t)=y(t)-C*xh_pc(:,t);
    P_pc{t+1}= A*P0_pc{t}*A'+V1;
    xf_pc(:,t)=xh_pc(:,t)+K0_pc{t}*e_pc(t);
    xh_pc(:,t+1) = A*xf_pc(:,t)+B*u(t);
end

norm(xh-xh_pc,inf)
norm(xf-xf_pc,inf)
