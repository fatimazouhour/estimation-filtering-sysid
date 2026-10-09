% WITHOUT GUI
clear all
close all
clc

load dryer2.mat

% load data make sure mean=0
u= u2-mean(u2);
y= y2-mean(y2);

u_m= mean(u);
y_m= mean(y);

% split up the data 
% we split 1.verification
%          2.estimation
Ntot= length(u);
Ne= Ntot/2;

ye= y(1:Ne);
yv= y(Ne+1:end);

ue= u(1:Ne);
uv= u(Ne+1:end);

% define datasets
ze= [ye,ue];
zv=[yv,uv];

% we can plot
figure
subplot(2,1,1), plot(1:Ntot,u,'b',1:Ne,ue,'g--',Ne+1:Ntot,uv,'r--'),ylabel('u'), grid on
subplot(2,1,2),plot(1:Ntot,y,'b',1:Ne,ye,'g--',Ne+1:Ntot, yv,'r--'),ylabel('y')
grid on, 

Nv= Ntot-Ne;
%--------------ARX MODEL----------------------
close all
for nk =1:4 
    for na=1:7 % usually not greater than 10
        nb=na;
        M= arx(ze,[na nb nk]); 
       
        %figure, resid(ze,M,'CORR',30);
        yh = compare(zv,M,1);
        %figure, plot(1:Nv,yv,'g',1:Nv,yh,'r')
       
        % validation criteria
        N0=10;  % usually first 7 but 10 is enough though
        MSE= 1/(Nv-N0)*norm(yv(N0+1:end)-yh(N0+1:end))^2;
        n=na+nb;
        n_arx(na)=n;
        FPE_arx(nk,na)= (Nv-N0+n)/(Nv-N0-n)*MSE;
        AIC_arx(nk,na) = n*2/(Nv-N0)+log(MSE);
        MDL_arx(nk,na)= n*log(Nv-N0)/(Nv-N0) + log(MSE);
        fit_arx(nk,na) =1-sqrt(MSE/(1/(Nv-N0)*norm(yv(N0+1:end)-mean(yv(N0+1:end)))^2));
        
       if nk ==3 & na==2, best_arx=M;
       end
    end
    
    %pause,close all
   
end
% -----notes on the previous
% nk=1, na=nb>=2 --> ok
% nk=2 , all OK
% nk=3, all OK
% nk=4, all OK 


% we can print & plot the performance
figure, plot(n_arx,FPE_arx(1,:),'r',n_arx,FPE_arx(2,:),'g',n_arx,FPE_arx(3,:),'b',n_arx,FPE_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('fpe arx')

figure, plot(n_arx,AIC_arx(1,:),'r',n_arx,AIC_arx(2,:),'g',n_arx,AIC_arx(3,:),'b',n_arx,AIC_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('aic arx')

figure, plot(n_arx,MDL_arx(1,:),'r',n_arx,MDL_arx(2,:),'g',n_arx,MDL_arx(3,:),'b',n_arx,MDL_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('mdl arx')

figure, plot(n_arx,fit_arx(1,:),'r',n_arx,fit_arx(2,:),'g',n_arx,fit_arx(3,:),'b',n_arx,fit_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('fit arx')


% to print
FPE_arx,AIC_arx,MDL_arx,fit_arx


present(best_arx)
[A,B,C,D,F] = polydata(best_arx);
A = tf(A,1,-1,'Variable','z^-1');
B = tf(B,1,-1,'Variable','z^-1'); 
TF_arx_pred_u = B
TF_arx_pred_y = 1-A

%--------ARMAX MODEL------------------------------
close all

for nk =1:4 
    for na=1:7
        nb=na;
        nc=na;
        M= armax(ze,[na nb nc nk]); 
       
        %figure, resid(ze,M,'CORR',30);
        yh = compare(zv,M,1);
        %figure, plot(1:Nv,yv,'g',1:Nv,yh,'r')
       
        % validation criteria
         N0=10;  % usually first 7 but 10 is enough though
         MSE= 1/(Nv-N0)*norm(yv(N0+1:end)-yh(N0+1:end))^2;
         n=na+nb+nc;
         n_armax(na)=n;
         FPE_armax(nk,na)= (Nv-N0+n)/(Nv-N0-n)*MSE;
         AIC_armax(nk,na) = n*2/(Nv-N0)+log(MSE);
         MDL_armax(nk,na)= n*log(Nv-N0)/(Nv-N0) + log(MSE);
         fit_armax(nk,na) =1-sqrt(MSE/(1/(Nv-N0)*norm(yv(N0+1:end)-mean(yv(N0+1:end)))^2));
         
    end
    
    %pause,close all
   
end


% notes about the model -->all ok 

figure, plot(n_armax,FPE_armax(1,:),'r',n_armax,FPE_armax(2,:),'g',n_armax,FPE_armax(3,:),'b',n_arx,FPE_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('fpe armax')

figure, plot(n_armax,AIC_armax(1,:),'r',n_armax,AIC_armax(2,:),'g',n_armax,AIC_armax(3,:),'b',n_arx,AIC_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('aic armax')

figure, plot(n_armax,MDL_armax(1,:),'r',n_armax,MDL_armax(2,:),'g',n_armax,MDL_armax(3,:),'b',n_arx,MDL_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('mdl armax')

figure, plot(n_armax,fit_armax(1,:),'r',n_armax,fit_armax(2,:),'g',n_armax,fit_armax(3,:),'b',n_arx,fit_arx(4,:),'k')
legend('nk=1','nk=2','nk=3','nk=4'), title('fit armax')


% to print
FPE_armax,AIC_armax,MDL_armax,fit_armax