%
% svd_notes.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Illustration of the geometrical meaning of the singular value
% decomposition, using a 2D example.
%

clc, clear, close all

bsvdgeometry = true;
bsvddiag = false;

if bsvdgeometry

% pick a matrix
%A = [1 1.5 ; 0 1]
A = randi([-2 2],2,2)

% singular value decomposition
[U,S,V] = svd(A)
U*S*V'

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

end

%--------------------------------------------------------------------------

if bsvddiag

% consider the effects of an svd on a diagonal matrix
% (e.g., the derivation in Aster, Eq 4.19, p. 100)
p = 6;
m = p+3;
n = p+2;
k = min(m,n);
s = 10*rand(p,1)
%G1 = diag(s);
G1 = diag(sort(s,'descend'));
G2 = [G1 ; zeros(m-p,p)];
G3 = [G1 zeros(p,n-p)];
G4 = [G1 zeros(p,n-p) ; zeros(m-p,n)];
G1, G2, G3, G4

% It is clear from here that G = S (as expected),
% and all the U and V are identiy matrices.
[U1,S1,V1] = svd(G1); whos G1 U1 S1 V1
[U2,S2,V2] = svd(G2); whos G2 U2 S2 V2
[U3,S3,V3] = svd(G3); whos G3 U3 S3 V3
[U4,S4,V4] = svd(G4); whos G4 U4 S4 V4

[Up1,Sp1,Vp1] = svd(G1,'econ'); whos G1 Up1 Sp1 Vp1, norm(G1 - Up1*Sp1*Vp1')
[Up2,Sp2,Vp2] = svd(G2,'econ'); whos G2 Up2 Sp2 Vp2, norm(G2 - Up2*Sp2*Vp2')
[Up3,Sp3,Vp3] = svd(G3,'econ'); whos G3 Up3 Sp3 Vp3, norm(G3 - Up3*Sp3*Vp3')

% I'm not sure why the econ command is not designed to work for Case 4.
% From "help svd": If X is m-by-n with m > n, then only the first
% n columns of U are computed and S is n-by-n.
[Up4,Sp4,Vp4] = svd(G4,'econ'); whos G4 Up4 Sp4 Vp4, norm(G4 - Up4*Sp4*Vp4')
% we can manually remove the last two columns
Sp4 = Sp4(1:p,1:p);
Up4(:,p+1:k) = [];
Vp4(:,p+1:k) = [];
whos G4 Up4 Sp4 Vp4, norm(G4 - Up4*Sp4*Vp4')

% note: you can take the inverse of these
whos Sp1 Sp2 Sp3 Sp4
inv(Sp1)
% note: you can't take the inverse of these
whos S1 S2 S3 S4

Gdag1 = Vp1*inv(Sp1)*Up1'
Gdag2 = Vp2*inv(Sp2)*Up2'
Gdag3 = Vp3*inv(Sp3)*Up3'
Gdag4 = Vp4*inv(Sp4)*Up4'

whos Gdag1 Gdag2 Gdag3 Gdag4
whos G1 G2 G3 G4

end

%--------------------------------------------------------------------------

