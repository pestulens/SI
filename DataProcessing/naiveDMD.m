close all;clc;clear
%% load data
DATA = readmatrix("20230908_40bar_60D_S84000.txt");
Data = DATA(11:end,:);
% Data = lowpass(DATA,5,1000);

t = Data(:,1) ;     %time (s)
Pt = Data(:,2) ;    %Tank pressure (bar)
Puv = Data(:,3);    %Pressure in control valve upstream (bar)
Pdv = Data(:,4);    %Pressure in control valve downstream (bar)
Pcat = Data(:,5);   %Pressure in catalyst (bar)
Pc = Data(:,6);     %Chamber pressure (bar)
T = Data(:,7);      %Temperature (degree)
dm = Data(:,8);     %Mass flow rate (g/s)
thru = Data(:,9);   %Thrust (kgf)
thru = thru-4.6;
x = find( Data(:,10)==1 ); % open valve signal
u = zeros(length(t),1);
u( x(1):x(end) ) = 100;
tu = x(1):x(end); %Epoch of open valve signal

% cal.t = t(x(1)-100:x(end));
% cal.Pt = Pt(x(1)-100:x(end));
% cal.Pcat = Pcat(x(1)-100:x(end));
% cal.PC = Pc(x(1)-100:x(end));
% cal.thru = abs(thru(x(1)-100:x(end)));
% cal.dm = dm(x(1)-100:x(end));
cal.t = t(x(1)+200:x(end));
cal.Pt = lowpass(Pt(x(1)+200:x(end)),1,1000);
% cal.Pt = 40+zeros(length(cal.t),1);
cal.Pcat = lowpass(Pcat(x(1)+200:x(end)),1,1000);
cal.PC = lowpass(Pc(x(1)+200:x(end)),1,1000);
cal.thru = abs(lowpass( thru(x(1)+200:x(end)),1,1000));
cal.dm = lowpass(dm(x(1)+200:x(end)),1,1000);

%% rearrange data
Data_r = [cal.Pt cal.Pcat cal.PC cal.thru]';
Data_r = Data_r(:,1:end-100);
% U = u( x(1)-100:x(end) );


%% DMD

% X   = [4 2 1 0.5;7 0.7 0.07 0.007];
% Xp  = [2 1 0.5 0.25;0.7 0.07 0.007 0.0007];
% Ups = [-4 -2 -1 -0.5] ;

X   = Data_r(:,1:end-1);
Xp  = Data_r(:,2:end);
approxA = X*pinv(Xp);

%% Simulation 
S.X(:,1) = Data_r(:,1);

% dt = 1/1000;
% for k=1:length(X)
%     S.X(:,k+1) = (approxA*S.X(:,k) + approxB*S.u(k));
% end
% Vf = @(t,x,u) approxA*x + approxB*u;

% [odet,odeX] = ode45(@(t,x) sys(t,x,approxA,approxB,Ups),0:0.001:5,Data_r(:,1));

tspan = 0:0.001:5;
x0 = [40;1;1;0];
[odet, odex] = ode45(@(odet,odex) sys(odet,odex,approxA),tspan,x0);

odex = odex';

%% plot 
% figure(1)
% plot(S.X(1,:)')
% title("Pressure in tank ")
% xlabel("epoch")
% ylabel("bar")
% 
% figure(2)
% plot(S.X(2,:)')
% title("Pressure in catalyst ")
% xlabel("epoch")
% ylabel("bar")
% 
% figure(3)
% plot(S.X(3,:)')
% title("Pressure in chamber")
% xlabel("epoch")
% ylabel("bar")
% 
% figure(4)
% plot(S.X(4,:)')
% title("Thrust")
% xlabel("epoch")
% ylabel("N")
% 
% figure(5)
% plot(S.X(5,:)')
% title("Mass flow rate")
% xlabel("epoch")
% ylabel("g/s")

figure(1)
plot(odex(1,:)')
title("Pressure in tank ")
xlabel("epoch")
ylabel("bar")

figure(2)
plot(odex(2,:)')
title("Pressure in catalyst ")
xlabel("epoch")
ylabel("bar")

figure(3)
plot(odex(3,:)')
title("Pressure in chamber")
xlabel("epoch")
ylabel("bar")

figure(4)
plot(odex(4,:)')
hold on
plot(Data_r(end,:))
title("Thrust")
xlabel("epoch")
ylabel("N")

figure(5)
plot(Ups)
title("Mass flow rate")
xlabel("epoch")
ylabel("g/s")


%% ode fun

function  dx = sys(t,x,A)
   dx = A*x ;
end
