clear all
close all
clc

load data1.mat

u1=u;
y1=y;

u = u-mean(u);
y= y-mean(y);

Ntot= length(u);
Ne= Ntot/2;
Nv= Ntot-Ne;

ue= u(1:Ne);
uv= u(Ne+1:end);

ye= y(1:Ne);
yv= y(Ne+1:end);

ze=[ye ue];
zv= [yv uv];

%% for ARX
for nk=1:4
    for na =1:7
        nb=na;
        
        model= arx(ze,[na nb nk]);
        %figure, resid(model,ze,'corr',30)
        yh= compare(zv,model,1);
        
        
        n=na+nb;
        n_arx(na)= n;
        N0=10;
        
        MSE= 1/(Nv-N0) * norm(yv(N0+1:end)-yh(N0+1:end))^2;
       
        RMSE_arx(nk,na)= sqrt(MSE);
        
        
        
    end
    %pause, close all
end
RMSE_arx
figure, plot(n_arx,RMSE_arx','o-'),title('RMSE_arx')
legend('nk=1','nk=2','nk=3','nk=4')
        

%% for ARMAX


for nk=1:4
    for na =1:7
        nb=na;
        nc=na;
        
        model= armax(ze,[na nb nc nk]);
     %   figure, resid(model,ze,'corr',30)
        yh= compare(zv,model,1);
        
        n=na+nb+nc;
        n_armax(na)= n;
        N0=10;
        
        MSE= 1/(Nv-N0) * norm(yv(N0+1:end)-yh(N0+1:end))^2;
       
        RMSE_armax(nk,na)= sqrt(MSE);
        
        if nk ==2& na==2, best_armax= model;  end     
    end
    %pause, close all
end

RMSE_armax
figure, plot(n_armax,RMSE_armax','o-'),title('RMSE_armax')
legend('nk=1','nk=2','nk=3','nk=4')

present(best_armax);
[A,B,C,D,F] = polydata(best_armax);
A= tf(A,1,-1,'variable','z^-1');
B= tf(B,1,-1,'variable','z^-1');
C= tf(C,1,-1,'variable','z^-1');

tf_u = B/C
tf_y= 1-A/C

%% for OE model 
for nk=1:4
    for nf =1:7
        nb=nf;
        
        
        model= oe(ze,[nf nb nk]);
       % figure, resid(model,ze,'corr',30)
        yh= compare(zv,model,1);
        
        n=nf+nb;
        n_oe(nf)= n;
        N0=10;
        
        MSE= 1/(Nv-N0) * norm(yv(N0+1:end)-yh(N0+1:end))^2;
       
        RMSE_oe(nk,nf)= sqrt(MSE);
        
        if nk==2 & nf==3 , best_oe = model;end
       
    end
  %  pause, close all
end

RMSE_oe 
figure, plot(n_oe,RMSE_oe','o-'),title('RMSE OE')
legend('nk=1','nk=2','nk=3','nk=4')

present(best_oe)
[A,B,C,D,F]= polydata(best_oe);

B= tf(B,1,-1,'variable','z^-1');
F= tf(F,1,-1,'variable','z^-1')

tf_u_oe= B/F;
%%
close all
clc

N = length(y1);

Phi = [-y(3:N-1) -y(2:N-2) -y(1:N-3) u(3:N-1) u(2:N-2) u(1:N-3)];
Y= y(4:N);
theta= Phi\Y

% EUI 2norm

epsilon2= 4 ; % given 
sigma= inv(Phi'*Phi);
theta_2min= theta-epsilon2*sqrt(diag(sigma));
theta_2max= theta+epsilon2*sqrt(diag(sigma));

EUI2= [theta_2min,theta_2max]


% to do the PUI2 we need to check that is it not empty 
alpha = norm(Y-Phi*theta);

if alpha^2 <=epsilon2^2
    theta_2min= theta-sqrt(diag(sigma))*sqrt(epsilon2^2-alpha^2);
     theta_2max= theta+sqrt(diag(sigma))*sqrt(epsilon2^2-alpha^2);
     PUI_2= [theta_2min, theta_2max];
end

alpha , PUI_2

epsilon_inf= 0.1;
A= inv(Phi'*Phi)*Phi';
for idx=1:length(theta),
    theta_min_inf(idx,1) = A(idx,:)*(Y-epsilon_inf*sign(A(idx,:))');
    theta_max_inf(idx,1) =    A(idx,:)*(Y+epsilon_inf*sign(A(idx,:))');
end

EUI_inf= [theta_min_inf,theta_max_inf]


%% 
close all
clear all
clc

