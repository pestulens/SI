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
cal.t = t(x(1)-100:x(end));
% cal.Pt = lowpass(Pt(x(1)-100:x(end)),5,1000);
cal.Pt = 40+zeros(length(cal.t),1);
cal.Pcat = lowpass(Pcat(x(1)-100:x(end)),5,1000);
cal.PC = lowpass(Pc(x(1)-100:x(end)),5,1000);
cal.thru = abs(lowpass( thru(x(1)-100:x(end)),5,1000));
cal.dm = lowpass(dm(x(1)-100:x(end)),5,1000);

%% rearrange data
Data_r = [cal.Pt cal.Pcat cal.PC cal.thru cal.dm]';
% U = u( x(1)-100:x(end) );


%% DMD
% 
% X   = [4 2 1 0.5;7 0.7 0.07 0.007];
% Xp  = [2 1 0.5 0.25;0.7 0.07 0.007 0.0007];
% Ups = [-4 -2 -1 -0.5] ;


X   = Data_r(:,1:end-1);
Xp  = Data_r(:,2:end);
Ups = u( x(1)-100:x(end)-1 )';

Omega = [X;Ups];

[U,Sig,V] = svd(Omega,'econ');

thresh = 1e-10;
rtil = length(find(diag(Sig)>thresh));

Util    = U(:,1:rtil); 
Sigtil  = Sig(1:rtil,1:rtil);
Vtil    = V(:,1:rtil); 

[U,Sig,V] = svd(Xp,'econ');

thresh = 1e-10;
r = length(find(diag(Sig)>thresh));

Uhat    = U(:,1:r); 
Sighat  = Sig(1:r,1:r);
Vbar    = V(:,1:r); 

n = size(X,1); 
q = size(Ups,1);
U_1 = Util(1:n,:);
U_2 = Util(n+q:n+q,:);

approxA = Uhat'*(Xp)*Vtil*inv(Sigtil)*U_1'*Uhat;
approxB = Uhat'*(Xp)*Vtil*inv(Sigtil)*U_2';


[W,D] = eig(approxA);

Phi = Xp * Vtil * inv(Sigtil) * U_1'*Uhat * W;

%% Simulation 
S.X(:,1) = Data_r(:,1);
S.u = Ups;
% dt = 1/1000;
% for k=1:length(X)
%     S.X(:,k+1) = (approxA*S.X(:,k) + approxB*S.u(k));
% end
% Vf = @(t,x,u) approxA*x + approxB*u;

% [odet,odeX] = ode45(@(t,x) sys(t,x,approxA,approxB,Ups),0:0.001:5,Data_r(:,1));
[odet, odex] = ode89(@(odet,odex) sys(odet,odex,approxA,approxB), 0:0.00001:5, [40;1;1;0.1;0.1]);
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
title("Thrust")
xlabel("epoch")
ylabel("N")

figure(5)
plot(odex(5,:)')
title("Mass flow rate")
xlabel("epoch")
ylabel("g/s")


%% ode fun

function  dx = sys(t,x,A,B)
   u = 100;
   dx = A*x + B*u;
end
