function [del_rng] = rngchn_mogi(n1,e1,depth,del_v,ning,eing,plook)
% Calculate a forward model for a Mogi source
%
%  USEAGE: [del_rng] = rngchn_mogi(n1,e1,depth,del_v,ning,eing,plook);
%
%  INPUT: n1 = local north coord of center of Mogi source (km)
%         e1 = local east coord of center of Mogi source (km)
%         depth = depth of Mogi source (km) for points to calculate
%                 range change. This vector is the depth of the 
%                 mogi source plus the elevation of the point
%                 taken from the DEM.
%         del_v = Volume change of Mogi source (km^3)
%         ning = north coord's of points to calculate range change
%         eing = east coord's of points to calculate range change
%
%  OUTPUT: del_rng = range change at coordinates given in ning and eing.
%                    If ning and eing are vectors with the same dimensions,
%                    del_rng is a vector. If ning is a row vector and eing
%                    is a column vecor, del_rng is a matrix of deformation
%                    values...compliant with Feigle and Dupre's useage.
%

[m,n] = size(ning);
[mm,nn] = size(eing);

%----coef for bulk modulus pressure <--> volume relation is below
%dsp_coef = (1000000*del_v*15)/(pi*16);

%----coef for McTigue's pressure <--> volume relation is below
dsp_coef = (1000000*del_v*3)/(pi*4);

if(mm==1 && n==1)
    %disp('Calculating a matrix of rngchg values')
    del_rng = zeros(m,nn);
    %del_d = del_rng;
    %del_f = del_rng;
    tmp_n = del_rng;
    tmp_e = del_rng;
    for i_loop = 1:m
        tmp_e(i_loop,:) = eing;
    end
    for i_loop = 1:nn
        tmp_n(:,i_loop) = ning;
    end
    % horizontal distance
    d_mat = sqrt((tmp_n-n1).^2 + (tmp_e-e1).^2);
    % denominator of displacement field for mogi source
    tmp_hyp = ((d_mat.^2 + depth.^2).^1.5);
    % Uh, horizontal displacement 
    del_d = dsp_coef*d_mat./tmp_hyp;
    % Uz
    del_f = dsp_coef*depth./tmp_hyp;
    % azimuthal angle
    azim = atan2((tmp_e-e1),(tmp_n-n1));
    % compute Ue and Un using Uh and azimuthal angle
    e_disp = sin(azim).*del_d;
    n_disp = cos(azim).*del_d;
    % project displacement field onto look vector
    for i_loop = 1:nn
        del_rng(:,i_loop) = ...
            [e_disp(:,i_loop) n_disp(:,i_loop) del_f(:,i_loop)]*plook';
    end
    % flip sign projected distance (or flip plook)
    del_rng = -1.0*del_rng;
    
elseif ((mm==1 && m==1) || (n==1 && nn==1))
    if (n~=nn)
        error('Coord vectors not equal length!')
    end
    %disp('Calculating a vector of rngchng values')
    % horizontal distance
    d_mat = sqrt((ning-n1).^2 + (eing-e1).^2);
    % denominator of displacement field for mogi source
    tmp_hyp = ((d_mat.^2 + depth.^2).^1.5);
    % Uh, horizontal displacement
    del_d = dsp_coef*d_mat./tmp_hyp;
    % Uz
    del_f = dsp_coef*depth./tmp_hyp;
    % azimuthal angle
    azim = atan2((eing-e1),(ning-n1));
    % compute Ue and Un using Uh and azimuthal angle
    e_disp = sin(azim).*del_d;
    n_disp = cos(azim).*del_d;
    % project displacement field onto look vector
    del_rng = [e_disp n_disp del_f]*plook';
    % flip sign of projected distance (or flip plook)
    del_rng = -1.0*del_rng;
    
else
    error('Coord vectors make no sense!')
end
