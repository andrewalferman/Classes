function [x,ProdData]= GenerateEnsembles2(n,M)
% n - grid size
% M - number of generated ensembles

data = load('hydrographdata.txt','ascii');
x = data(:,1);  
y = data(:,2);
a = y(2);
t0 = x(1);
t1 = x(end);

% figure
% plot(data(:,1),data(:,2),'--b','LineWidth',2); hold on;
% pause

% ft = fittype(sprintf('%0.7f+(b*c/(sqrt(pi)*(x-s)))*exp(-c^2*(log(abs(x-s))-d)^2)',a),...
%     'coefficients',{'s','b','c','d'},...
%     'independent','x');
% 
% f = fit(x,y,ft,'Startpoint',[12 1000 1.5 4]);
% 
% ConfidenceLevel=.99; % 99% confidence interval
% ConfInt = (confint(f,ConfidenceLevel))';

load('confint'); % provides ConfInt from file
ConfInt(3,:)=(2)*ConfInt(3,:); % change correlation length here
ConfInt(1,:)=(3)*ConfInt(1,:); % add variance in time
x = linspace(t0,t1,n);
ProdData = zeros(length(x),M);
rng('default');

for j=1:M
U = rand(4,1);
r = U.*(ConfInt(:,2)-ConfInt(:,1))+ConfInt(:,1); %generate coefficients for the model
ProdData(:,j) = a+(r(2)*r(3))./(sqrt(pi)*(x-r(1))).*exp(-r(3)^2*(log(abs(x-r(1)-50))-r(4)).^2);
end

