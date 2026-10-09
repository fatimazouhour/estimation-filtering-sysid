clear all
close all
clc

load heater.mat

ue= ue-mean(ue);
uv=uv-mean(uv);

ye=ye-mean(ye);
yv=yv-mean(yv);

Ne= length(ue);
Nv= length(uv);

ze= [ye,ue];
zv= [yv,uv];


%----ARX----------
for nk=1:5
    for na=1:5
        nb=na;
        
        model= arx(ze,[na,nb,nk]);
    %    figure, resid(ze,model,'corr',30)
    yh= compare(zv,model,1);
    
    n= na+nb;
    n_arx(na)= n;
    
    N0=10;
    MSE= 1/(Nv-N0)*norm(yv(N0+1:end)-yh(N0+1:end))^2;
    
    MDL_arx(nk,na)= n*log(Nv-N0)/(Nv-N0)+log(MSE);
    
    RMSE_arx(nk,na)= sqrt(MSE);
    if nk==4 & na==2, 
        bestmodel= model;
    end
    
    end
   % pause, close all
end

MDL_arx,RMSE_arx

figure, plot(n_arx,MDL_arx','o-'),title('mdl ARX')
legend('nk=1','nk=2','nk=3','nk=4','nk=5')
figure, plot(n_arx,RMSE_arx','o-'),title('RMSE ARX')
legend('nk=1','nk=2','nk=3','nk=4','nk=5')

%------------ARMAX-----------

for nk=1:5
    for na=1:5
        nb=na;
        nc= na;
        
       model= armax(ze,[na,nb,nc,nk]);
     % figure, resid(ze,model,'corr',30)
     yh= compare(zv,model,1);
     n= na+nb+nc;
    n_armax(na)= n;
    
    N0=10;
    MSE= 1/(Nv-N0)*norm(yv(N0+1:end)-yh(N0+1:end))^2;
    
    MDL_armax(nk,na)= n*log(Nv-N0)/(Nv-N0)+log(MSE);
    
    RMSE_armax(nk,na)= sqrt(MSE);
    end
 %   pause, close all
end
RMSE_armax,MDL_armax


figure, plot(n_armax,MDL_armax','o-'),title('mdl ARMAX')
legend('nk=1','nk=2','nk=3','nk=4','nk=5')
figure, plot(n_armax,RMSE_armax','o-'),title('RMSE ARMAX')
legend('nk=1','nk=2','nk=3','nk=4','nk=5')

%----------OE------------
for nk=1:5
    for nf=1:5
        nb=nf;
        
        model= oe(ze,[nf,nb,nk]);
    %    figure, resid(ze,model,'corr',30)
        
    end
   % pause, close all
end


present(bestmodel),
[A,B,C,D,F]= polydata(bestmodel);
A=tf(A,1,-1,'variable','z^-1')
B= tf(B,1,-1,'variable','z^-1')

F2= 1-A
F1=B




   