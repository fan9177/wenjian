% Modal analysis CD-head

% Mass matrix
M1 = 3e-4; % [kg] focus coil
M2 = 5e-4; % [kg] radial coil
M3 = 2e-4; % [kg] lens

M = [
    M1 0  0;
    0  M2 0;
    0  0  M3];

% Stiffness matrix
freq1 = 40; % [Hz] eigenfreq lensassy;

K02 = (2*pi*freq1)^2 * (M2+M1+M3);
K12 = 0.2e6; % [N/m]
K23 = 1.0e6; % [N/m]

K = [
      K12 -K12           0;
     -K12 (K02+K12+K23) -K23;
      0   -K23           K23];

% Determine eigenvectors, eigenvalues 
[V,D] = eig(K,M);
D = diag(D);
n = size(V);

% % Scale to max norm to bring the mass matrix to 'realistic' values
% for i=1:n(2)
%     V(:,i) = V(:,i)/max(V(:,i));
% end

% Eigenvalues to eigenfrequencies
D = sqrt(D)/(2*pi);

% Reduced modal mass and stiffness matrices

MM = V'*M*V;
KK = V'*K*V;

% Take the diagonal elements

MM = diag(MM);
KK = diag(KK);

% Modal damping
zeta = 0.01;
CC = 2*zeta*sqrt(MM.*KK);

% Transfer from Force point A to Displacement point B
A = 2;
B = 3;
for i=1:n(2)
    mode(i) = tf([V(A,i)*V(B,i)],[MM(i) CC(i) KK(i)]);
end

% Plot response
figure(1); bode( mode(1), mode(2), mode(3), logspace(1,5,1000)*2*pi);
figure(2); bode( mode(1), 'r', mode(2), 'r', mode(3), 'r', mode(1)+mode(2)+mode(3), 'k', logspace(1,5,1000)*2*pi);
