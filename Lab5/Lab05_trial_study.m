clear all
close all
clc

a1= 0.93;
b1= 1.5;
b2=-3;

theta0= [a1 b1 b2]';

rng('default')

Nmax=500;
u = 10*2*rand(Nmax,1)-10;


for t= 1:Nmax
    e(t)= 3*randn;
    if t==1
        phi(:,t) =[0 0 0]';
    elseif t==2
        phi(:,t)= [-y(t-1) u(t-1) 0]';
    else 
        phi(:,t)=[-y(t-1) u(t-1) u(t-2)]';
    end
    y(t,1) = phi(:,t)'*theta0+e(t);
    
end

figure,plot(y), title('y')

Phi= [];
Rtol = 1e-10;

for t=1:Nmax
    if t==1
        phi(:,t) =[0 0 0]';
    elseif t==2
        phi(:,t)= [-y(t-1) u(t-1) 0]';
    else 
        phi(:,t)=[-y(t-1) u(t-1) u(t-2)]';
    end
    Phi=[Phi;phi(:,t)'];
    Y= y(1:t);
    Rt= Phi'*Phi/t;
    if abs(det(Rt))>Rtol
        theta(:,t)= Phi\Y;
    else 
        theta(:,t)=zeros(size(theta0));
    end
end

theta_n = theta(:,Nmax)
theta_n , theta0


figure, hold on
plot(theta(1,:),'r','linewidth',1.5)
plot(theta(2,:),'g','linewidth',1.5)
plot(theta(3,:),'k','linewidth',1.5)
for idx=1:length(theta0)
    refline(0,theta0(idx))
end

grid  on, title('est params for arx(1,2,1)')
legend('a1','b1','b2')

% andersons whiteness test 
for t= 1:Nmax
    if t==1
        phi(:,t) =[0 0 0]';
    elseif t==2
        phi(:,t)= [-y(t-1) u(t-1) 0]';
    else 
        phi(:,t)=[-y(t-1) u(t-1) u(t-2)]';
    end
    yh(t,1) = phi(:,t)'*theta_n; % with the one we estimated 
    
    
end
e_pred= y-yh;
e_pred_mean= mean(e_pred);

figure,plot(1:Nmax,y,'b',1:Nmax,yh,'r--') , grid on, xlim([0 100])
legend('output','pred y')

[r,lags]= xcorr(e_pred,'coeff');
figure, plot(lags,r,'ro-'),xlim([0,50])
refline(0,2/sqrt(Nmax)),refline(0,-2/sqrt(Nmax))



