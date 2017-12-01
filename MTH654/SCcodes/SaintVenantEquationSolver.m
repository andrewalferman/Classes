function [TimeSteps,YChosen,QChosen,WSEChosen,VChosen] = SaintVenantEquationSolver(tdata,Qdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a,b] - space domain
% Nx - number of grid points in the space domain including the endpoints
% tdata in [0,T] - time domain
% Qdata - inflow hydrograph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options=optimset('Display','off');

%% data setup
Nt = length(tdata); % number of time steps
% time steps for illustration purposes only, actual time step should satisfy CFL condition
T=tdata(end); % assume uniform time step, i.e., tdata = 0:T/Nt:T; 
dt=tdata(2)-tdata(1);

%figure; plot(tdata, Qdata, '--bo');
%hold on

% Extend time interval to see flow througout channel
T=2*T;

tt=[tdata,tdata(end)+dt:dt:T]; 
Qdatarep=[Qdata;repmat(Qdata(end),length(tt)-Nt,1)];% repeat last value to extrapolate

% Build spline to interpolate at arbitrary t
Qspline=spline(tt, Qdatarep); 


%% domain setup
a = 0; % left endpoint of the domain
ChannelLength = 500; % length of the reach
t0 = 0; % initial time

%% numerical scheme parameters
Nx = 20; % number of grid points in the space domain
Cn = 1; % Courant number
theta = 0.8; % weighting coefficient in the numerical scheme
dx = ChannelLength/(Nx-1); % grid size
x = (a:dx:(a+ChannelLength))'; % grid points

%% channel parameters setup
% ChannelWidth (B) - width of the channel
% ConvFactor (k) - conversion factor:  (Length^(1/3)/Time), 1 m^(1/3)/s for SI, or 1.4859 ft^(1/3)/s U.S. customary units, if required. 
% GMCoef (n) - Gauckler-Manning coefficient (defines the channel material) 
% Reasonable values are given at http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm

ChannelWidth = 50; % width of the channel
B=ChannelWidth;
GMCoef = 0.035; % Gauckler-Manning coefficient
ConvFactor = 1; % conversion factor
g = 9.80665; % gravity constant
S0 = 0.0005; % channel slope

%% setup initial conditions

%initialize (arbitrary)
Q0=ppval(Qspline,t0);
%WSEn = WSE(t0)*ones(Nx,1); % water surface elevation
%Yn = WSEn - S0*(a+ChannelLength-x); % water depth
%Y0=Yn(1);
% Vn = Q0/(ChannelWidth*Yn(1,1))*ones(Nx,1); % flow velocity
Qn = Q0*ones(Nx,1); % flow discharge constant

% solve for steady state
    %Y' = (Sb-Sf) / (1 - (Q/A)^2/g/Y ) & A=BY & Q=VA & R=A/P & P=2Y+B & Sf=
    %n^2*Q^2 / k^2/A^2/R^(4/3) .
    syms y(t); A=B*y; P=2*y+B; R=A/P; Sf=GMCoef^2*Q0^2/ConvFactor^2/A^2/R^(4/3);
    f=(S0-Sf) / (1-(Q0/A)^2/g/y);
    odef = odeFunction(f,y);
    [x,Ys] = ode45(odef,x,9.7518);

    %A = Q0^2*GMCoef^2;
    %C = Q0^2/(ChannelWidth^2*g);
    %[x,Ys] = ode45(@(x,Y) (S0 - (A*(ChannelWidth+2*Y)^(4/3))/(ChannelWidth^(10/3)*Y^(10/3)))/(1-C/Y^3), x, Y0);

% use steady state as initial condition    
Yn=Ys;
Yn(end)=10; % Enforce boundary condition exactly
Vn = Q0./(ChannelWidth*Yn); % flow velocity
WSEn=Yn+S0*(a+ChannelLength-x);

%% calculation of Froude number
Fr = max((Qn./(abs(Yn).^(3/2)))/(ChannelWidth*sqrt(g)));
% disp(['Froude number is ',num2str(Fr),'.']);
% if Fr < 1
%     disp('Flow is subcritical.');
% elseif Fr > 1
%     disp('Flow is supercritical.');
% else
%     disp('Flow is critical.');
% end


%% additional outputs
location = [Nx]; % location of outflows
%location = [ceil(Nx/2),Nx]; % location of outflows
%location = 1:1:length(x);
l = length(location); % number of chosen locations
VChosen(1,1:l) = Vn(location)';
YChosen(1,1:l) = Yn(location)';
WSEChosen(1,1:l) = (Yn(location)+S0*(a+ChannelLength-x(location)))';
TimeSteps(1) = t0;

%% plot initial conditions for water stage, discharge and velocity
figsave=gcf;
% figure(9);
% subplot(3,1,1); plot(x,WSEn); 
% title('Initial water stage');
% xlabel('x');
% ylabel('WSE');
% 
% subplot(3,1,2); plot(x,Qn); 
% title('Initial discharge'); 
% xlabel('x');
% ylabel('Discharge');
% 
% subplot(3,1,3); plot(x,Vn); 
% title('Initial velocity'); 
% xlabel('x');
% ylabel('Velocity');
% pause(0.5);

%% constants of the discretization scheme
c(1) = theta/dx;
c(2) = (1-theta)/dx;
c(3) = ChannelWidth; 
c(4) = GMCoef; % n
c(5) = ConvFactor; % k
c(6) = theta;
c(7) = 1-theta;
c(8) = (c(4)/c(5))^2;
c(9) = g;

%% cycle setup
t = t0;
nt = 1; % time step index
UGuess = [Yn;Vn]; % inital guess for the first time step

%% initial time step
celerity = sqrt(abs(g*Yn)); 
dt = 0.97*Cn*min(dx./(Vn + celerity)); % time step should satisfy CFL condition
% disp(['Initial time step is ',num2str(dt)]);

%% cycle over all time steps
while t < T
    
    t = t + dt;
    
    DiffYn = Yn(2:1:Nx,1) - Yn(1:1:Nx-1,1);
    SumYn = Yn(1:1:Nx-1,1) + Yn(2:1:Nx,1);
    
    DiffVn = Vn(2:1:Nx,1) - Vn(1:1:Nx-1,1);
    SumVn = Vn(1:1:Nx-1,1) + Vn(2:1:Nx,1);
    
    Sfn = c(8)*4^(1/3)*abs(SumVn).*SumVn./((c(3)*SumYn./(SumYn + c(3))).^(4/3));
    
    NonLinearSystem = @(U) DiscreteSystem(U,t,DiffYn,SumYn,DiffVn,SumVn,S0,Sfn,c,dt,Qspline);
    
    [U,fval,exitflag,output]  = fsolve(NonLinearSystem,UGuess,options);
%     disp('The value of the function:'); disp(fval);    
%     disp(['Exit flag is ', num2str(exitflag)]);
    
    if exitflag ~= 1
        
        t = t - dt;
        dt = dt/2;
        
        if nt ~= 1
            UGuess = dt*(Un - Unminus1)/dtn + Un;
        end  
        
    else
        
        Unminus1 = [Yn;Vn];
        Un = U;
        
        dtn = dt;
%         disp(['Time step is ',num2str(dt)]);
        Yn = U(1:Nx,1);
        Vn = U(Nx+1:2*Nx,1);
       
        celerity = sqrt(abs(g*Yn)); 
        dt = abs(0.97*Cn*min(dx./(Vn + celerity))); % time step should satisfy CFL condition
        
        UGuess = dt*(Un - Unminus1)/dtn + Un;
        
        nt = nt + 1;
        TimeSteps(nt) = t;
        VChosen(nt,1:l) = Vn(location)';
        YChosen(nt,1:l) = Yn(location)';
        WSEChosen(nt,1:l) = YChosen(nt,1:l) + S0*(a+ChannelLength-x(location))';

        WSEn = Yn + S0*(a+ChannelLength-x); % water surface elevation
        Qn = ChannelWidth*Vn.*Yn;

        

%         figure(10);
%         subplot(3,1,1); plot(x,WSEn); 
%         title(['Water stage at time t = ' num2str(t) ', dt = ' num2str(dt)]);
%         xlabel('x');
%         ylabel('WSE');
%         
%         subplot(3,1,2); plot(x,Qn); 
%         title(['Discharge at time t = ' num2str(t) ', dt = ' num2str(dt)]); 
%         xlabel('x');
%         ylabel('Discharge');
%         
%         subplot(3,1,3); plot(x,Vn); 
%         title(['Velocity at time t = ' num2str(t) ', dt = ' num2str(dt)]); 
%         xlabel('x');
%         ylabel('Velocity');
%         
%         pause(0.5)
    end
end

%% output the depth and flow at chosen location
figure(figsave)
subplot(3,1,1); plot(TimeSteps',WSEChosen);hold on 
title('Water stage at chosen location'); 
%legend('250 m');
xlabel('Time'); 
ylabel('WSE');
% subplot(3,1,1); plot(TimeSteps',YChosen);hold on 
% title('Water height at chosen location'); 
% %legend('250 m');
% xlabel('Time'); 
% ylabel('Y');

QChosen = ChannelWidth*VChosen.*YChosen;
subplot(3,1,2); plot(TimeSteps',QChosen); hold on 
title('Discharge at chosen location'); 
%legend('250 m', '500 m');
xlabel('Time'); 
ylabel('Discharge Q');

subplot(3,1,3); plot(TimeSteps',VChosen); hold on 
title('Velocity at chosen location'); 
%legend('250 m', '500 m');
xlabel('Time'); 
ylabel('Velocity V');

pause(0.1) % force figure to refresh

end


function F = DiscreteSystem(U,t,DiffYn,SumYn,DiffVn,SumVn,S0,Sfn,c,dt,option)
    
    Nx = length(U)/2;
    
    Y = U(1:Nx,1);
    V = U(Nx+1:2*Nx,1);
    
    F = ones(2*Nx,1);
    
    DiffY = Y(2:1:Nx,1) - Y(1:1:Nx-1,1);
    SumY = Y(1:1:Nx-1,1) + Y(2:1:Nx,1);
    
    DiffV = V(2:1:Nx,1) - V(1:1:Nx-1,1);
    SumV = V(1:1:Nx-1,1) + V(2:1:Nx,1);
    
    AveY = 0.5*(c(6)*SumY + c(7)*SumYn);
    AveV = 0.5*(c(6)*SumV + c(7)*SumVn);
    
    dYdx = c(1)*DiffY + c(2)*DiffYn;
    dVdx = c(1)*DiffV + c(2)*DiffVn;
    
    Sf = c(8)*4^(1/3)*abs(SumV).*SumV./((c(3)*SumY./(SumY + c(3))).^(4/3));
    AveSf = 0.5*(c(6)*Sf + c(7)*Sfn);
    
    F1 = SumY - SumYn + 2*dt*(AveY.*dVdx + AveV.*dYdx);
    F2 = SumV - SumVn + 2*dt*(AveV.*dVdx + c(9)*dYdx + c(9)*(AveSf-S0));
    
    F(1,1) = c(3)*V(1,1)*Y(1,1) - Discharge(t,option); % VBy = Q at t=t^(n+1)
    F(2:2:2*Nx-2,1) = F1;
    F(3:2:2*Nx-1,1) = F2;
    F(2*Nx,1) = Y(Nx,1) - Depth(t);
end

function y = Depth(t)
    y = 10;
end

function y = WSE(t)
    y = 10;
end

function Q = Discharge(t,Qspline)
Q=ppval(Qspline,t);

end
        
function d = Depth_a(y,YOld,QOld,Qa,ChannelWidth)

    Q = QOld(2,1);
    Y = YOld(2,1);
    
    d = Qa/(ChannelWidth*y)-2*sqrt(abs(9.8*y))-(Q/(ChannelWidth*Y)-2*sqrt(abs(9.8*Y)));
end
