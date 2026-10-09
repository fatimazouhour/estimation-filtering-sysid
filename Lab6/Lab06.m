clear 
clc 

load heater.mat

% the data sets are already split 

% yet we have to remove the mean value

ue_m = ue -mean(ue) ;
uv_m = uv- mean(uv);

ye_m= ye-mean(ye);
yv_m= yv-mean(yv);

% split the data

ze= [ye_m, ue_m];
zv= [yv_m, uv_m];

% 

Ne= length(ue);
Nv= length(uv);


% then start with ARX
for nk= 1:5
    for na=1:5
        nb=na;
        model= arx(ze,[na, nb,nk]);
        
   % figure, resid(ze,model,'corr',30);
    y= compare(zv,model,1); % since for prediction
    N0= 10;
    %figure, plot(1:Nv,yv_m,'-k',1:Nv,y,'--r')
    MSE_arx= 1/(Nv-N0)*norm(yv_m(N0+1:end)-y(N0+1:end))^2; 
    RMSE_arx(nk,na)= sqrt(MSE_arx);
    
    n(na)= na+nb; % for the arx model 
    MDL_arx(nk,na)= n(na)* log(Nv-N0)/(Nv-N0)+log(MSE_arx);
   
    end
  % pause , close all
end

%-----------notes on graph-----------


MDL_arx, RMSE_arx

% figure,
% plot(n, MDL_arx,'o-'),title('mdl')
% legend('nk1','nk=2','nk=3','nk=4','nk=5')
% 
% figure,
% plot(n,RMSE_arx,'o-'), title('rmse')
% legend('nk1','nk=2','nk=3','nk=4','nk=5')



%----------OE model------------------
for nk= 1:5,
    for nb= 1:5
        nf=nb;
        
        model=oe(ze,[nb nf nk]);
        
        %figure, resid(ze, model,'corr',30);
        y= compare(zv,model,1);
        
        n(nb)= nb+nf;
        
         MSE_oe= 1/(Nv-N0)*norm(yv_m(N0+1:end)-y(N0+1:end))^2; 
         RMSE_oe(nk,nb)= sqrt(MSE_oe);
         MDL_oe(nk,nb)= n(nb)* log(Nv-N0)/(Nv-N0)+log(MSE_oe);
    end
    %pause, close all
end


MDL_oe, RMSE_oe


% plot(n, MDL_oe,'o-'),title('mdl')
% legend('nk1','nk=2','nk=3','nk=4','nk=5')
%  
% figure, plot(n,RMSE_oe,'o-'), title('rmse')
% legend('nk1','nk=2','nk=3','nk=4','nk=5')
% WE CAN TELL THAT OE IS NOTTT A GOOD CANDIDATE FROMTHE START

%------------------NARX----------------------

Netdet = ['HHHHHH';'L-----'];
na= 1;
nb= na;
nk=1;

trparms= settrain;
trparms= settrain(trparms,'maxiter',500);
[W1,W2]= nnarx(Netdet,[na,nb,nk],[],[],trparms,ye_m',ue_m');
figure, drawnet(W1,W2)
    
alpha= W2(1:end-1);
alpha_0=W2(end) ;

beta=W1(:,1:end-1);
beta_0 = W1(:,end);


% to get the predicted output 
tmin= max([na+1,nb+nk]);

for t=tmin:Nv
    phi= [yv_m(t-1:-1:t-na);uv_m(t-nk:-1:t-nk-nb+1)];
    yp(t)=alpha*tanh(beta*phi+beta_0)+alpha_0;
end

figure, plot(1:Nv,yv_m,'g',1:Nv,yp,'r')

MSE_narx= 1/(Nv-N0)*norm(yv_m(N0+1:end)-yp(N0+1:end)')^2; 

RMSE_narx(nk,na)= sqrt(MSE_narx);
n= na+nb;
MDL_narx(nk,na)= n(na)* log(Nv-N0)/(Nv-N0)+log(MSE_arx);
