clear
close all
clc

rng('default');

Nmax=500; % given
u= rand(Nmax,1)*10*2-10; % see explanation on paper

% ARX parameters : coefficients
a1=0.93;
b1= 1.5;
b2=-3;

theta_0=[a1,b1,b2]';

% ARX output computation

for t=1:Nmax
   e(t) = 3*randn;
   if t==1
       phi(:,t) = [0 0 0]';   
   elseif t==2
       phi(:,t)=[-y(t-1) u(t-1) 0]';
   else
       phi(:,t)=[-y(t-1) u(t-1) u(t-2)]';
   end
   
   y(t,1) = phi(:,t)'*theta_0+e(t);
end

% plot the ARX output
figure, plot(y) , title('system output')
           
% ------------part2--------------------------------

% see from slide 52
Phi= [];% initialize empty matrix
R_tol= 1e-10;

for t=1:Nmax
    if t==1
       phi(:,t) = [0 0 0]';   
   elseif t==2
       phi(:,t)=[-y(t-1) u(t-1) 0]';
   else
       phi(:,t)=[-y(t-1) u(t-1) u(t-2)]';
   end
    Phi= [Phi;phi(:,t)']; %in every iteration a new vector of phi is added to the big matrix
    Y=y(1:t); % initialize Y matrix
    R_t= Phi'*Phi*1/t;
    
    if abs(det(R_t))>R_tol
        theta_t(:,t)= Phi\Y; % note fina n3ml compare m3 method 8er 
    else 
        theta_t(:,t)=zeros(size(theta_0));
    end
end

theta_N= theta_t(:,Nmax)

% we can plot these estimates
figure, hold on
plot(theta_t(1,:))
plot(theta_t(2,:))
plot(theta_t(3,:))

%-------------Part3-----------------------------------------

for t= 1:Nmax
    if t==1
       phi(:,t) = [0 0 0]';   
    elseif t==2
       phi(:,t)=[-y(t-1) u(t-1) 0]';
    else
       phi(:,t)=[-y(t-1) u(t-1) u(t-2)]';
    end
    y_h(t,1) = phi(:,t)'*theta_N;
end

eps= y-y_h; % prediction error
eps_mean= mean(eps);

figure,
plot(1:Nmax,y,'b',1:Nmax,y_h,'r--'),grid on 
xlim([0,100])
legend('output','predicted output')

figure,
plot(eps,'r'),grid on, title('eps')
refline(0,eps_mean)

% anderson's
[r,lags]= xcorr(eps,'coeff'); % memorize!!
figure,
plot(lags,r,'ro-'),xlim([0,50])
refline(0,2/sqrt(Nmax))
refline(0,-2/sqrt(Nmax))
    

