%SCdriver -- driver file for applying Stochastic Collocation to Saint
%Venant System with KL inflows

%% Get KL expansion that retains 99% of variance
[t, MeanK, lambda, phi, myNKL]=KLinflow(.99);
% Returns t, MeanK, lambda, phi, myNKL (number of terms to keep)

%% Find Stochastic Collocation nodes (assuming uniform [-sqrt(3),sqrt(3)])

% Choose Clenshaw-Curtis 1D level:
level=3;

% Full Tensor Clenshaw-Curtis:
[xf,wf]=cc_grid_dataset(myNKL,2^(level-1)+1); % defined by number of points per dimension

% Sparse Grid Clenshaw-Curtis:
[xs,ws]=sparse_grid_cc_dataset(myNKL,level-1); % level off by one since they start at zero

% plot first two dimensions of grid:
figure(3); clf;
plot(xf(1,:),xf(2,:),'o')
hold on
plot(xs(1,:),xs(2,:),'o','MarkerFaceColor','blue')
legend('Full','Sparse')
title(['Clenshaw-Curtis grids, L= ', num2str(level)])
pause(0.1)

%% Perform Discrete Projection:

% we will want simulations upto 2*T
tt=[t(1):t(2)-t(1):2*t(end)];


%% Full Tensor:
tic
MaxQf=0; % initialize expected max of Qout
MeanQf=0*tt; % initialize expected Qout
MeanQf2=0*tt; % intialize expected Qout^2
MaxTf=0; % initialize expected time of max Qout

figure(4);clf;
figure(5);clf;
for i=1:size(xf,2) 
    % use parfor here (if you have it) for speedup, but note that plots
    % within loop won't display
  
   % Build Qin by evaluating KL at collocation point xf(:,i)
   Qin=MeanK; 
   for j=1:myNKL
       Qin=Qin+sqrt(3)*xf(j,i)*sqrt(lambda(j))*phi(:,j);
   end
    figure(4);plot(t,Qin);hold on
    
   % Simulate system at collocation point
   figure(5);
   [t2,Y,Qout,WSE,V] = SaintVenantEquationSolver(t,Qin);
   
   % A Quantity of Interest is peak outflow
   [maxQ,imax]=max(Qout);
   MaxQf = MaxQf + wf(i)*maxQ; % add according to weight

   % Need to find time of max outflow for problem 3
   MaxTf = MaxTf + wf(i)*t2(imax); %add according to weight
   
   % mean & std
   Qspline=spline(t2,Qout);
   MeanQf=MeanQf+wf(i)*ppval(Qspline,tt);
   MeanQf2=MeanQf2+wf(i)*ppval(Qspline,tt).^2;

   
   pause(0.1) % force figures to refresh
end
figure(4);title('Inflow collocation points full grid')


disp(['Expected MaxQf=',num2str(MaxQf)])
disp(['Expected MaxTf=',num2str(MaxTf)])
toc



%% Repeat for Sparse:
tic
MaxQs=0;
MaxTs=0;
figure(6);clf;
figure(7);clf;
for i=1:size(xs,2)
  
   % Build Qin by evaluating KL at collocation point xs(:,i)
   Qin=MeanK;

   for j=1:myNKL
       Qin=Qin+sqrt(3)*xs(j,i)*sqrt(lambda(j))*phi(:,j);
   end
   figure(6);plot(t,Qin);hold on
   
   % Simulate system at collocation point
   figure(7)
   [t2,Y,Qout,WSE,V] = SaintVenantEquationSolver(t,Qin);
   
   % A Quantity of Interest is peak outflow
   [maxQ,imax]=max(Qout);
   MaxQs = MaxQs + ws(i)*maxQ; % add according to weight

   % Need to find time of max outflow for problem 3
   MaxTs = MaxTs + ws(i)*t2(imax); %add according to weight
   
   pause(0.1) % force figures to refresh
end
figure(6);title('Inflow collocation points sparse grid')

disp(['Expected MaxQs=',num2str(MaxQs)])
disp(['Expected MaxTf=',num2str(MaxTs)])
toc

StdDev = 0*tt;
StdPlus=0*tt;
StdMinus=0*tt;
StdDev=ppval(Qspline,tt)-MeanQf;
StdPlus=ppval(Qspline,tt)+2*StdDev;
StdMinus=ppval(Qspline,tt)-2*StdDev;

figure(8);
plot(t2,Qout);
hold on;
plot(tt,StdPlus,'--');
plot(tt,StdMinus,'--');
legend('Expected Outflow', 'Plus 2 Std Dev.', 'Minus 2 Std. Dev.');
%plot(tt,MeanQf2);
title('Expected Outflow with Standard Deviations');

