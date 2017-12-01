function [t, MeanK, lambda, phi, myNKL]=KLinflow(energytol)
n = 200; % time points
M = 10;  % number of data ensembles
NKL = 5;  % maxiumum number of terms to keep in the KL

%% Input the data (although here we are generating synthetic data)
[t,ProdData] = GenerateEnsembles2(n,M);
ProdData=ProdData+100; %shift to be realistic

figure(1)
plot(t,ProdData,'b');
title(['Synthesized Data, M= ',num2str(M)]);
pause(0.1)

%% Compute Mean and Covariance of data
MeanK = mean(ProdData,2);
C=cov(ProdData')';

%% Compute the eigenvalues and eigenfunctions of the data
[lambda,phi] = IntEqSolver1d(t(1),t(end),length(t),C,NKL,'collocation',1,0);

%% Compute energies
energyF=zeros(NKL,1);
sumlambda=sum(lambda);
for i=1:NKL
energyF(i)=sum(lambda(1:i))/sumlambda;
end

%% Display results
% figure;
% plot(lambda,'-ok');
% title('Eigenvalues');
% 
% figure;
% plot(t,phi,'Linewidth',2);
% LegendBase=['1st';'2nd';'3rd';'4th';'5th';'6th'];
% Legend={};
% for i=1:NKL
%     Legend{i}=strcat(LegendBase(i,:),' (',num2str(100*energyF(i)),'%)');
% end
% legend(Legend);   
% title('Eigenfunctions');

nv=find(energyF>energytol); % energytol percent of energy retained
myNKL=nv(1);

figure(2); clf
plot(t,MeanK,'--k','Linewidth',4);
hold on
plot(t,phi(:,1:myNKL)*diag(sqrt(lambda(1:myNKL))),'Linewidth',2);
legend('Mean','1st','2nd','3rd','4th','5th','6th');
title(['Scaled Eigenfunctions, NKL= ',num2str(myNKL)]);
pause(0.1)

%% sample the KL to generate simulated realizations
if (0)
    
    nrel = 20;
    dK = zeros(n,nrel);
    for j = 1:nrel
        % Z = sqrt(lambda).*randn(NKL,1); % assume random coefficients are N(0,1)
        % dK(:,j) = phi*Z;
        
        Z = sqrt(lambda(1:myNKL)).*randn(myNKL,1); % assume random coefficients are N(0,1)
        dK(:,j) = phi(:,1:myNKL)*Z;
        
    end
    
    figure(2)
    plot(t,dK+MeanK*ones(1,size(dK,2)));
    hold on
    hm=plot(t,MeanK,'--k','Linewidth',4);
    legend(hm,'Mean');
    title(['Sampled realizations with myNKL=',num2str(myNKL)]);
end