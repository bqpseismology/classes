%
% svd_notes.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Misc concepts associated with the singular value decomposition
%

clc, clear, close all

bsvdgeometry = true;    % geometrical meaning of SVD for 2D example
bsvddiag = false;       % svd for a diagonal matrix

%--------------------------------------------------------------------------

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

%print(gcf,'-depsc','svd_2D');

end

%--------------------------------------------------------------------------

if bsvddiag

% consider the effects of an svd on a diagonal matrix
% (e.g., the derivation in Aster, Eq 4.19, p. 100)
% [This is a fairly technical point, but some of the commands below may be
% illuminating in general.]

p = 6;      % fixed for all cases
m = p+3;
n = p+2;
k = min(m,n);
s = 10*rand(p,1)
%D1 = diag(s);
% consider four cases from Aster, Ch. 3
D1 = diag(sort(s,'descend'));
D2 = [D1 ; zeros(m-p,p)];
D3 = [D1 zeros(p,n-p)];
D4 = [D1 zeros(p,n-p) ; zeros(m-p,n)];
D1, D2, D3, D4
disp(sprintf('D1 is %i (m) x %i (n) with p = %i',size(D1),p));
disp(sprintf('D2 is %i (m) x %i (n) with p = %i',size(D2),p));
disp(sprintf('D3 is %i (m) x %i (n) with p = %i',size(D3),p));
disp(sprintf('D4 is %i (m) x %i (n) with p = %i',size(D4),p));

% It is clear from here that D = S (as expected),
% and all the U and V are identity matrices.
% note: why does it switch the order of the last two basis vectors in U4 and V4?
[U1,S1,V1] = svd(D1); %whos D1 U1 S1 V1
[U2,S2,V2] = svd(D2); %whos D2 U2 S2 V2
[U3,S3,V3] = svd(D3); %whos D3 U3 S3 V3
[U4,S4,V4] = svd(D4); %whos D4 U4 S4 V4

% economy size svd
[Up1,Sp1,Vp1] = svd(D1,'econ'); %whos D1 Up1 Sp1 Vp1, norm(D1 - Up1*Sp1*Vp1')
[Up2,Sp2,Vp2] = svd(D2,'econ'); %whos D2 Up2 Sp2 Vp2, norm(D2 - Up2*Sp2*Vp2')
[Up3,Sp3,Vp3] = svd(D3,'econ'); %whos D3 Up3 Sp3 Vp3, norm(D3 - Up3*Sp3*Vp3')

% I'm not sure why the econ command is not designed to work for Case 4.
% From "help svd": If X is m-by-n with m > n, then only the first
% n columns of U are computed and S is n-by-n.
[Up4,Sp4,Vp4] = svd(D4,'econ'); %whos D4 Up4 Sp4 Vp4, norm(D4 - Up4*Sp4*Vp4')
% we can manually remove the last two columns
Sp4 = Sp4(1:p,1:p);
Up4(:,p+1:k) = [];
Vp4(:,p+1:k) = [];
%disp('manually corrected version:');
%whos D4 Up4 Sp4 Vp4, norm(D4 - Up4*Sp4*Vp4')

% note: you can take the inverse of these
whos Sp1 Sp2 Sp3 Sp4
disp('inverse of Sp:');
inv(Sp1)
% note: you can't take the inverse of these
whos S1 S2 S3 S4

Ddag1 = Vp1*inv(Sp1)*Up1'
Ddag2 = Vp2*inv(Sp2)*Up2'
Ddag3 = Vp3*inv(Sp3)*Up3'
Ddag4 = Vp4*inv(Sp4)*Up4'

whos Ddag1 Ddag2 Ddag3 Ddag4 D1 D2 D3 D4

whos Up1 Up2 Up3 Up4 Vp1 Vp2 Vp3 Vp4

end

%--------------------------------------------------------------------------

