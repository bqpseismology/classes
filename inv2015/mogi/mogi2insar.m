function rad_calc_rot2 = mogi2insar(x,y,z,v,iplot,imask)
%MOGI2INSAR generate a synthetic interferogram based on a Mogi source
%
% Mogi source is located at [x, y, z] with a volume change of v;
% the calculated LOS range changes (in mm) are stored in rad_calc_rot2
% if plot_ind==1, then plot the simulated interferogram 
%

%global eing_vec ning_vec obs_phavec plook

bvc = [x y z v 0 0 0 0]';

% INPUT TRACK AND LOOK ANGLES 
% satellite track angle from the geographic north (in degrees) 
track =  -13.3;
% radar look angle from the vertical (in degrees) 
look  = 23.0;

% INPUT FILE WIDTH (SAMPLE), LENGTH (LINE), AND POSTING
SAMPLE = 1100;
LINE = 980;
POSTING = 40;

%===========================

track = track*pi/180.0;
look = look*pi/180.0;
% radar line-of-sight (LOS) vector in [E, N, U]
% note: there is a sign flip in rngchn_mogi.m that could instead be applied here
plook = [-sin(look)*cos(track) sin(look)*sin(track) cos(look)];

num_rows = LINE;
num_cols = SAMPLE;

m = num_rows;
n = num_cols;

%---Below is temp. for something Zhong wants
TMP = ones(m,n);

%---posting is POSTING
ning = 0:POSTING:(m-1)*POSTING;
eing = 0:POSTING:(n-1)*POSTING;

%-----change m to km
ning = ning/1000;
eing = eing/1000;

[eing_vec_bg,ning_vec_bg,vec_ind_bg] = get_ind_pha(TMP,ning,eing);
calc_mat = zeros(size(TMP));

%bvc = parm_srt(1:7,1)

% calculate forward model: LOS dispalcement in mm
calc_rng = rngchn_mogi(bvc(2),bvc(1),bvc(3),bvc(4),ning_vec_bg,eing_vec_bg,plook);
%calc_mat(vec_ind_bg) = (calc_rng - eing_vec_bg*bvc(6) - ning_vec_bg*bvc(7) - bvc(5));
calc_mat(vec_ind_bg) = calc_rng;

%rng2rad = (2.0*pi)/28.3;
rng2rad = 1;        % (do not apply)
rad_calc = calc_mat*rng2rad;
rad_calc_rot = (flipud(rad_calc))';
% pha_file = [filename,'.syn'];
% fid = fopen(pha_file,'w', 'b');
% fwrite(fid,rad_calc_rot,'float');
% fclose(fid);

% whether you want to plot with the data mask
if imask==1
    load data_mask      % saved from read_data.m
    rad_calc_rot2 = rad_calc_rot'.*data_mask;
    rad_calc_rot2(rad_calc_rot2==0) = nan;
else
    rad_calc_rot2 = rad_calc_rot';
end

if iplot==1
    plot_model(rad_calc_rot2);
end
