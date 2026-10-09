clear all
close all
clc

load dryer2.mat

u=u2-mean(u2);
y=y2-mean(y2);

N= length(u);
Ne= N/2;
Nv= N-Ne;

ue= u(1:Ne);
uv= u(Ne+1:end);

ye= y(1:Ne);
yv= y(Ne+1:end);

ze= [ye,ue];
zv=[yv,uv];

%------------ARX--------------------------
for nk=1:4
    for na= 1:7
        nb= na;
        
        model = arx(ze,[na,nb,nk]);
    %   figure, resid(ze,model,'corr',30)
        yh= compare(zv,model,1);
         
        n= na+nb;
        n_arx(na)= n;
        
        N0=10;
        MSE= 1/(Nv-N0)*norm(yv(N0+1:end)-yh(N0+1:end))^2;
        
        RMSE_arx(nk,na)= sqrt(MSE);
        FPE_arx(nk,na)= (Nv-N0+n)/(Nv-N0-n)*MSE;
        AIC_arx(nk,na)= n*2/(Nv-N0)+log(MSE);
        MDL_arx(nk,na)= n*log(Nv-N0)/(Nv-N0)+log(MSE);
        
        bestfit_arx(nk,na)= 1-sqrt(MSE/(1/(Nv-N0)))*1/norm(yv(N0+1:end)-mean(y(N0+1:end)));
        
        
          if nk ==3 & na==2 , 
        bestmodel= model;
    end
    end
    
    RMSE_arx,FPE_arx,AIC_arx,MDL_arx,bestfit_arx
     
  
    
    
   % pause,close all
end

figure, plot(n_arx,RMSE_arx','o-'),title('RMSE arx')
legend('nk=1','nk=2','nk=3','nk=4')

figure, plot(n_arx,FPE_arx','o-'),title('FPE arx')
legend('nk=1','nk=2','nk=3','nk=4')


figure, plot(n_arx,AIC_arx','o-'),title('AIC arx')
legend('nk=1','nk=2','nk=3','nk=4')

figure, plot(n_arx,MDL_arx','o-'),title('MDL arx')
legend('nk=1','nk=2','nk=3','nk=4')

figure, plot(n_arx,bestfit_arx','o-'),title('best fit arx')
legend('nk=1','nk=2','nk=3','nk=4')
%------------ARMAX-----------------------
for nk=1:4
    for na= 1:7
        nb= na;
        nc= na;
        
        model = armax(ze,[na,nb,nc,nk]);
   %     figure, resid(ze,model,'corr',30)
   
        yh= compare(zv,model,1); 
         n= na+nb+nc;
        n_armax(na)= n;
        
        N0=10;
        MSE= 1/(Nv-N0)*norm(yv(N0+1:end)-yh(N0+1:end))^2;
        
        RMSE_armax(nk,na)= sqrt(MSE);
        FPE_armax(nk,na)= (Nv-N0+n)/(Nv-N0-n)*MSE;
        AIC_armax(nk,na)= n*2/(Nv-N0)+log(MSE);
        MDL_armax(nk,na)= n*log(Nv-N0)/(Nv-N0)+log(MSE);
        
        bestfit_armax(nk,na)= 1-sqrt(MSE/(1/(Nv-N0)))*1/norm(yv(N0+1:end)-mean(yv(N0+1:end)));
        
    end
  %  pause,close all
end

RMSE_armax,FPE_armax,AIC_armax,MDL_armax,bestfit_armax

figure, plot(n_armax,AIC_armax','o-'),title('AIC armax')
legend('nk=1','nk=2','nk=3','nk=4')

figure, plot(n_armax,MDL_armax','o-'),title('MDL armax')
legend('nk=1','nk=2','nk=3','nk=4')

figure, plot(n_armax,bestfit_armax','o-'),title('best fit armax')
legend('nk=1','nk=2','nk=3','nk=4')

%--------------OE---------------------------
for nk=1:4
    for nf= 1:7
        nb= nf;
        
        model = oe(ze,[nf,nb,nk]);
   %    figure, resid(ze,model,'corr',30)
        
    end
  % pause,close all
end

%%
        
present(bestmodel),
[A,B,C,D,F]= polydata(bestmodel);
A= tf(A,1,-1,'variable','z^-1')
B= tf(B,1,-1,'variable','z^-1')

tf_u  = B
tf_y= 1-A


