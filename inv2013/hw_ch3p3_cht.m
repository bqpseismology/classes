%
% hw_ch3p3_cht.m
% Problem 3.3
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
%
% exploration of SVD spaces and solutions using parabola fits to sets of points
% 

clc
clear
close all

% true parameters (m1, m2, m3) in the form of Example 1.1
m_true = [10 20 9.8]';

% forward model (in-line function)
G = @(t) ([ones(size(t)) t -0.5*t.^2]);

% finely discretized time vector for plotting solutions
tmin = 0.9;
tmax = 3.1;
t_vec = linspace(tmin,tmax,100);
d_vec = @(m) ( m(1) + m(2)*t_vec - 0.5*m(3)*t_vec.^2);

% number of solutions to plot
avec = -5:5;
na = length(avec);

%---------

disp('======== (a) ========');
% (a) forward problem matrix for full, consistent data set, three data points
t1 = [1 2 3]';
G1 = G(t1);
y_true = G1*m_true;
y = y_true;             % noise free predicted data
[U1,S1,V1] = svd(G1);

% solution
p = rank(G1);
%m_recov1 = inv(G1)*y
%m_recov1 = pinv(G1)*y
m_recov1 = (V1*inv(S1)*U1')*y

figure; subplot(2,2,1);
plot(t1,y,'o',t_vec,d_vec(m_recov1),'r');
xlabel('Time (s)'); ylabel('Height (m)'); xlim([tmin tmax]);
title({'(a) Three data points:','a unique exact solution'})

disp(['p = ' num2str(p)]);

%---------

disp('======== (b) ========');
% (b) nonunique case, two data points (nontrivial model null space)
t2 = [1 2]';
G2 = G(t2);
y_true = G2*m_true;
y = y_true;          % noise free predicted data
[U2,S2,V2] = svd(G2);

% solution
p = rank(G2);
Vp = V2(:,1:p);
Sp = S2(1:p,1:p);
Up = U2(:,1:p);
m_recov2 = Vp*inv(Sp)*Up'*y
%m_recov2 = pinv(G2)*y

subplot(2,2,2); hold on;
% add in scaled null space vectors here to show nonuniqueness
for ii=1:na
    alpha = avec(ii);
    m = m_recov2 + alpha*V2(:,3);
    plot(t2,y,'o',t_vec,d_vec(m),'r');
    mnorm(ii) = norm(m);
    resid(ii) = norm(G2*m-y);
end
xlabel('Time (s)'); ylabel('Height (m)'); xlim([tmin tmax]);
title({'(b) Two data points:','infinite number of exact solutions'})

disp(['p = ' num2str(p)]);
disp('model null space vector:')
null_m = V2(:,3)

%---------

disp('======== (c) ========');
% (c) inconsistent case (nontrivial data null space)
t3 = [1 2 2.1 3]';
G3 = G(t3);
y_true = G3*m_true;
% predicted data (with noise)
y = y_true;
y(2) = y_true(2)+1;
y(3) = y_true(3)-1;
[U3,S3,V3] = svd(G3);

% solution
p = rank(G3);
Vp = V3(:,1:p);
Sp = S3(1:p,1:p);
Up = U3(:,1:p);
m_recov3 = Vp*inv(Sp)*Up'*y
%m_recov3 = pinv(G3)*y

subplot(2,2,3);
plot(t3,y,'o',t_vec,d_vec(m_recov3),'r');
xlabel('Time (s)'); ylabel('Height (m)'); xlim([tmin tmax]);
title({'(c) Four data points:','one least-squares (approximate) solution'})
disp(['p = ',num2str(p)])
disp('data null space vector:')
null_d = U3(:,4)

%---------

disp('======== (d) ========');
% (d) both null spaces nontrivial
t4 = [2 2]';     % TWO MEASUREMENTS AT SAME TIME
G4 = G(t4);
y_true = G4*m_true;
% predicted data (with noise)
y = y_true;
y(1) = y_true(1)+1;
y(2) = y_true(2)-1;
[U4,S4,V4] = svd(G4);

% solution
p = rank(G4);
Vp = V4(:,1:p);
Sp = S4(1:p,1:p);
Up = U4(:,1:p);
m_recov4 = Vp*inv(Sp)*Up'*y
%m_recov4 = pinv(G4)*y

subplot(2,2,4); hold on;
% add in scaled null space vectors here to produce a suite of equally good models,
% and this show nonuniqueness of the least-squares solution.
for ii=1:na
    alpha = avec(ii);
    m = m_recov4 + alpha*(10*V4(:,2)+V4(:,3));
    plot(t4,y,'o',t_vec,d_vec(m),'r');
    mnorm(ii) = norm(m);
    resid(ii) = norm(G2*m-y);
end
xlabel('Time (s)'); ylabel('Height (m)'); xlim([tmin tmax]);
title({'(d) Two data points: infinite number of',
    'least-squares (approximate) solutions'})

% display p and the model space null vectors
disp(['p = ' num2str(p)])
disp('data null space vector:')
null_d = U3(:,4)
disp('model null space vectors:')
null_m_1 = V3(:,2)
null_m_2 = V3(:,3)

%if iprint==1, orient tall; print(gcf,'-dpsc',sprintf('%shw_ch1',pdir)); end

%==========================================================================
