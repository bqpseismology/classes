%
% test function for computing 2D norms
%

function [f,ind] = testfun2D(x,y,p)

[m,n] = size(x);
f = zeros(m,n);
imat = zeros(m,n);

fin = zeros(m*n,1);
xin = reshape(x,m*n,1);
yin = reshape(y,m*n,1);

% find which domain each input point belongs to
i1 = find( and(xin >= 8, yin >= 8) );
i2 = find( and(xin >= 12, yin <= 4) );
i3 = find( and(xin <= 2, yin <= 2) );
i4 = find( and(xin <= 1, yin >= 15) );

% everything NOT in the four squares is in this domain
i5 = setdiff([1:m*n],unique([i1 ; i2; i3; i4]));

% save indexing
imat(i1) = 1; imat(i2) = 2; imat(i3) = 3; imat(i4) = 4; imat(i5) = 5;

% assign function values
fin(i1) = p(1); fin(i2) = p(2); fin(i3) = p(3); fin(i4) = p(4); fin(i5) = p(5);

f = reshape(fin,m,n);
ind = reshape(imat,m,n);

%------------------

if 0==1
    pvec = [-2 -6 8 4 1]';
    % multiple of 16 OR nx --> infty will give the same value
    nxvec = [16 32 50 100 150 200 400];
    nf = length(nxvec);
    nvec = zeros(nf,1);
    
    for ii = 1:nf
        nx = nxvec(ii);
        xvec = linspace(0,16,nx);
        yvec = linspace(0,16,nx);
        [X,Y] = meshgrid(xvec,yvec);
        n = nx*nx;
        dx = xvec(2) - xvec(1);

        [F,Ind] = testfun2D(X,Y,pvec);
        f = F(:);

        if 0==1
            fn(1) = length( find(f == -2 ) ) / n;
            fn(2) = length( find(f == -6 ) ) / n;
            fn(3) = length( find(f == 3 ) ) / n;
            fn(4) = length( find(f == 4 ) ) / n;
            fn(5) = length( find(f == 1 ) ) / n;
            fn'

            figure; hold on;
            %pcolor(X,Y,F); shading flat;
            imagesc(F);
            %plot(X(:),Y(:),'.')
            axis equal, axis([0 nx+1 0 nx+1])
            colorbar
        end

        %Atot = n*dx^2;               % NOTE: does not equal xran*yran
        Atot = 16^2;
        dA = Atot/n;
        dAvec = ones(n,1)*dA;
        Avec = sqrt( dAvec/Atot );
        n*dA - Atot, sum(Avec.^2)         % check

        nvec(ii) = sum( f.^2 .* Avec.^2 );
    end
    
    % so the norm is about 5 -- now consider two custom meshes
    nx = 100;
    xvec = linspace(0,16,nx);
    yvec = linspace(0,16,nx);
    [X,Y] = meshgrid(xvec,yvec);
    n = nx*nx;
    dx = xvec(2) - xvec(1);
    [F,Ind] = testfun2D(X,Y,pvec);
    f = F(:);
    icen = (Ind==5);
    xcen = X(icen);
    ycen = Y(icen);
    fcen = F(icen);
    figure; plot(xcen,ycen,'.');
    
    % assuming a regular mesh
    %Atot = n*dx^2;
    Atot = 16^2;
    dA = Atot/n;
    dAvec = ones(n,1)*dA;
    Avec = sqrt( dAvec/Atot );
    n*dA - Atot, sum(Avec.^2)         % check
    sum( f.^2 .* Avec.^2 )
    
    % now add the four custom points with the correct associated area
    % NOTE: it doesn't matter WHERE they are for the norm calculation --
    %       just so that the area associated with each point is correct
    dAvec1 = dAvec(icen);
    dAvec1 = [dAvec1 ; 8^2 ; 4^2 ; 2^2 ; 1^2];
    Avec1 = sqrt( dAvec1/Atot );
    f1 = [fcen ; pvec(1:end-1) ];
    sum(Avec1.^2)
    sum( f1.^2 .* Avec1.^2 )
    
    % now a mesh with five "points"
    A5 = 16^2 - sum([8^2 ; 4^2 ; 2^2 ; 1^2]);
    dAvec2 = [8^2 ; 4^2 ; 2^2 ; 1^2 ; A5];
    Avec2 = sqrt( dAvec2/Atot );
    f2 = pvec;
    sum(Avec2.^2)
    sum( f2.^2 .* Avec2.^2 )
end
    
%----------------------------------------
