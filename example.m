clear all;
close all;
clc;
%% Original system
% A = [1.5 0; 0 0.1];
% B = [1;0];
% X(:,1) = [4;7];
% u = [-4 -2 -1 -0.5 -0.25 -0.125 -0.0625];
% for i =1:length(u)
%     X(:,i+1) = A*X(:,i)+B.*u(i);
% end

%% System Identification

X   = [4 2 1 0.5 0.25 0.125 0.0625;7 0.7 0.07 0.007 0.0007 0.00007 0.000007];
Xp  = [ 2 1 0.5 0.25 0.125 0.0625 0.03125;0.7 0.07 0.007 0.0007 0.00007 0.000007 0.0000007];
Ups = [-4 -2 -1 -0.5 -0.25 -0.125 -0.0625];
X_t = [4 2 1 0.5 0.25 0.125 0.0625 0.03125;7 0.7 0.07 0.007 0.0007 0.00007 0.000007 0.0000007];


% 
% X   = StateData(:,1:end-1);
% Xp  = StateData(:,2:end);
% Ups = InputData(:,1:end-1);

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

%% system reform

S.X(:,1) = [4;7];
S.Y = Ups;
for i =1:length(X)
    S.X(:,i+1) = approxA*S.X(:,i)+approxB.*S.Y(i);
end
Rerror(1,:) = abs(X_t(1,:)-S.X(1,:))./X_t(1,:);
Rerror(2,:) = abs(X_t(2,:)-S.X(2,:))./X_t(2,:);

Aerror(1,:) = abs(X_t(1,:)-S.X(1,:));
Aerror(2,:) = abs(X_t(2,:)-S.X(2,:));

%% plot 
figure(1)
plot(X_t(1,:))
hold on
plot(S.X(1,:))
hold on 
legend("original 1","SI 1")

figure(2)
plot(X_t(2,:))
hold on
plot(S.X(2,:))
legend("original 2","SI 2")
% legend("original 1","original 2","SI 1","SI 2")
