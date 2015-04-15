%
% lab_pca_matlab.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Principal Component Analysis of diets in Europe.
%
% TWO KEY OPERATIONS TO THE COLUMNS OF THE DATA MATRIX
% (1) centered: subtract mean
% (2) standardized (or scaled): divide by the standard deviation
%
% use the function pca; princomp is depricated and should not be used
%

clc
clear
close all
format compact

idata = 2;
bmatlabfigs = false;

%-------------------------------------------------
% from matlab
% http://www.mathworks.com/help/stats/feature-transformation.html#f75476

switch idata
    case 1
        load cities, whos, X = ratings;
    case 2
        X = load('~/GEOTOOLS/classes/inv2015/data/protein_matlab.dat');
end
[n,p] = size(X);

%-------------------------------------------------
% from matlab
% http://www.mathworks.com/help/stats/feature-transformation.html#f75476

if and(bmatlabfigs,idata==1)

categories
figure()
boxplot(ratings,'orientation','horizontal','labels',categories)
C = corr(ratings,ratings);
w = 1./var(ratings);
[wcoeff,score,latent,tsquared,explained] = pca(ratings,'VariableWeights',w);
coefforth = inv(diag(std(ratings)))*wcoeff;
cscores = zscore(ratings)*coefforth;

% visualizing the results
figure()
plot(score(:,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

figure()
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

[st2,index] = sort(tsquared,'descend'); % sort in descending order
extreme = index(1);
names(extreme,:)

biplot(coefforth(:,1:2),'scores',score(:,1:2),'varlabels',categories);
axis([-.26 0.6 -.51 .51]);

%figure()
%biplot(coefforth(:,1:3),'scores',score(:,1:3),'obslabels',names);
%axis([-.26 0.8 -.51 .51 -.61 .81]);
%view([30 40]);

end

%-------------------------------------------------
% from wiki: http://en.wikipedia.org/wiki/Principal_component_analysis

% REMOVE THE MEAN FOR EACH COLUMN
u = mean(X);
h = ones(n,1);
% centered matrix
B = X - h*u;
%B = X - repmat(u,n,1);

% COVARIANCE MATRIX
C = (1/(n-1)) * (B') * B;
%norm(C - cov(B))       % check
%norm(C - cov(X))       % check

% standard deviations of each column of X
% note: column vector diag(C) is equivalent to row vector var(X)
s = sqrt(diag(C))'
%s = sqrt(var(X))

Cdiag = diag(diag(C));
hCdiag = sqrtm(Cdiag);

% CORRELATION MATRIX (is a scaled covariance matrix)
%R = corrcoef(X);
%R = corrcoef(B);
%R = corrcov(C);
%R = C./(s'*s);
%R = diag(1./s)*C*diag(1./s);
R = inv(hCdiag)*C*inv(hCdiag)
Ccheck = hCdiag*R*hCdiag;
norm(Ccheck - C)

% STANDARDIZED (AND CENTERED) MATRIX
%Z = B ./ (h*s)
%Z = B * diag(1./s)
Z = B*inv(hCdiag);
norm(Z - zscore(X))
% check reverse operation
%Bcheck = Z .* (h*s');
%norm(B - Bcheck)

% EIGENVALUE DECOMPOSITION
[Vc,Dc] = eig(C);
%norm(C*Vc - Vc*Dc)

% if needed: sort eigenvalues and rearrange V
eigval = diag(Dc);
[~,isort] = sort(eigval,'descend');
Vc = Vc(:,isort);
Dc = diag(eigval(isort));
%norm(C*Vc - Vc*Dc)
eigval = diag(Dc);

[Vr,Dr] = eig(R);
% Vr * Dr * Vr'
% inv(hCdiag)*Vc * Dc * (inv(hCdiag)*Vc)'

if 0==1
    %%
    A = rand(8,3)
    s = sqrt(var(A))
    h = ones(8,1);
    % this looks like dividing each column by its standard deviation
    A ./ (h*s)
    A ./ repmat(s,8,1)
    % this is the series of elemetary matrix operations
    A * diag([1/s(1) 1 1]) * diag([1 1/s(2) 1]) * diag([1 1 1/s(3)])
    % this is what it means
    A * diag(1./s)
end

disp('mean(X), mean(B), mean(Z):');
mean(X), mean(B), mean(Z)
disp('var(X), var(B), var(Z):');
var(X), var(B), var(Z)
disp('std(X), std(B), std(Z):');
std(X), std(B), std(Z)

% SINGULAR VALUE DECOMPOSITION of B
% the scores matrix is nothing more than U*S from the SVD
% Vb = Vc (allowing sign changes)
[Ub,Sb,Vb] = svd(B);
svalb = diag(Sb);
USb = Ub*Sb
Vb
% check
norm( svalb.^2/(n-1) - eigval )

% SINGULAR VALUE DECOMPOSITION of Z
% Vz = Vr (allowing sign changes)
[Uz,Sz,Vz] = svd(Z)
svalz = diag(Sz);
USz = Uz*Sz
Vz

% Test 1: use centered matrix
% VB = Vb = Vc (allowing for some sign flips on columns of V)
% USB = Ub*Sb
[VB,USB] = pca(B)
% this is equivalent, since pca will center the matrix (i.e., remove mean)
%[V,US] = pca(X)
Bcheck = USB * VB';
norm(B - Bcheck)
% orthonormal:
norm(VB'*VB - eye(p))
Bcheck = USB * VB';
norm(B - Bcheck)
US1_check = B*VB;
norm(USB - US1_check)

% Test 2 (example used in matlab)
% note: this gives different US and V from Test 1
w = 1./(s.^2);
[Vw,USw] = pca(B,'VariableWeights',w);
% this is equivalent
%[Vw,USw] = pca(B,'VariableWeights','variance')
Bcheck = USw * Vw';
norm(B - Bcheck)
disp('Vw is NOT orthornomal:');
norm(Vw'*Vw - eye(p))
Vworth = inv(hCdiag)*Vw;   % note: inv(hCdiag) = diag(1./s)
% orthonormal:
norm(Vworth'*Vworth - eye(p))
% this shows how USw can be computed
USw_check = Z*Vworth;
norm(USw - USw_check)

% Test 3: use centered+standardized matrix as input
% USZ = USw (allowing for sign flips)
% VZ = Vz = Vr = inv(hCdiag)*Vw  (allowing for sign flips)
[VZ,USZ] = pca(Z);
Zcheck = USZ * VZ';
norm(Z - Zcheck)
Bcheck = USZ * VZ' * hCdiag;
norm(B - Bcheck)
% orthonormal:
norm(VZ'*VZ - eye(p))
%VZ_check = inv(hCdiag)*Vw

%==========================================================================
