%
% svd_notes.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Illustration of the geometrical meaning of the singular value
% decomposition, using a 2D example.
%

clc, clear, close all

% pick a matrix
A = [1 1.5 ; 0 1]
%A = randi([-2 2],2,2)

% singular value decomposition
[U,S,V] = svd(A)

s1 = S(1,1);
s2 = S(2,2);
v1 = V(:,1);
v2 = V(:,2);
u1 = U(:,1);
u2 = U(:,2);

% unit circle for plotting
R = 1;
n = 100;
t = linspace(0,2*pi,n);
x = R*cos(t);
y = R*sin(t);
Vx = [x ; y];

ax0 = 3*[-1 1 -1 1];
figure; nr=1; nc=2; f = 1.2; fsize = 12;

subplot(nr,nc,1); hold on;
plot(Vx(1,:),Vx(2,:),'b-'); axis equal, axis(ax0);
plot([0 v1(1)],[0 v1(2)],'b',[0 v2(1)],[0 v2(2)],'b','linewidth',2);
text(f*v1(1),f*v1(2),'v_1','fontsize',fsize,'horizontalalignment','center');
text(f*v2(1),f*v2(2),'v_2','fontsize',fsize,'horizontalalignment','center');

Ux = A*Vx;      % for plotting
su1 = s1*u1;
su2 = s2*u2;

subplot(nr,nc,2); hold on;
plot(Ux(1,:),Ux(2,:),'r'); axis equal, axis(ax0);
%plot([0 u1(1)],[0 u1(2)],'r',[0 u2(1)],[0 u2(2)],'r','linewidth',2);
%text(f*u1(1),f*u1(2),'u_1','fontsize',fsize,'horizontalalignment','center');
%text(f*u2(1),f*u2(2),'u_2','fontsize',fsize,'horizontalalignment','center');
plot([0 su1(1)],[0 su1(2)],'r',[0 su2(1)],[0 su2(2)],'r','linewidth',2);
text(f*su1(1),f*su1(2),'\sigma_1 u_1','fontsize',fsize,'horizontalalignment','center');
text(f*su2(1),f*su2(2),'\sigma_2 u_2','fontsize',fsize,'horizontalalignment','center');