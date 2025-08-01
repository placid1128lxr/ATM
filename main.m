clc
clear
close all

global rot_e 
global v_light
global dt_data
global dataset
global mask_angle SNR_lim
global init_lat init_lon init_h

global ARAIM_URA_GPS ARAIM_URA_GLO ARAIM_URA_GAL ARAIM_URA_BDU ...
       ARAIM_BIAS_GPS ARAIM_BIAS_GLO ARAIM_BIAS_GAL ARAIM_BIAS_BDU

global ARAIM_URE_GPS ARAIM_URE_GLO ARAIM_URE_GAL ARAIM_URE_BDU ...
       ARAIM_BIAS_CONT_GPS ARAIM_BIAS_CONT_GLO ARAIM_BIAS_CONT_GAL ARAIM_BIAS_CONT_BDU;

global ARAIM_PSAT_GPS ARAIM_PSAT_GAL ARAIM_PSAT_GLO ARAIM_PSAT_BDU

global ARAIM_SIN_USRMASK_GPS ARAIM_SIN_USRMASK_GLO ARAIM_SIN_USRMASK_GAL ...
       ARAIM_SIN_USRMASK_BDU

global SIG_ACC_MAX_VERT
global SIG_ACC_MAX_HOR1 SIG_ACC_MAX_HOR2
global VPLT HPLT EMTT ATTEMPT_OPT

global ARAIM_USRMASK_GPS ARAIM_USRMASK_GLO ARAIM_USRMASK_GAL ARAIM_USRMASK_BDU 

global ARAIM_PCONST_GPS ARAIM_PCONST_GAL ARAIM_PCONST_GLO ARAIM_PCONST_BDU 

global PHMI_VERT PHMI_HOR P_THRES PFA_VERT PFA_HOR P_EMT PL_TOL FC_THRES ...
    PL0_FDE FDE_FLAG FDE_WF_FLAG

global flag_sig2_if

global T_Global flag_opt

global GDOP

rot_e = 7.2921151467e-5;
v_light = 299792458;
micro_g_to_meters_per_second_squared = 9.80665E-6;

dt_data = 1;
mask_angle = 15;
SNR_lim = 15;
flag_sig2_if = 0;
FDE_FLAG = 0;

T_Global=chi2inv(1-1e-7,1:30)'; 

GNSS_config.epoch_interval = 1;
GNSS_config.intend_no_GNSS_meas = 30;
GNSS_config.mask_angle = mask_angle;
GNSS_config.mask_SignalStrenth = SNR_lim;
GNSS_config.ISBBDS_GPS = 12.1*v_light/1e9;
GNSS_config.omit = [];
GNSS_config.GPSPRNList=1:37;
GNSS_config.BDSPRNList=174:210;

TC_KF_config.init_att_unc = deg2rad(20);
TC_KF_config.init_vel_unc = 0.1;
TC_KF_config.init_pos_unc = 5;
TC_KF_config.init_b_a_unc = 10000 * micro_g_to_meters_per_second_squared;
TC_KF_config.init_b_g_unc = deg2rad(10) / 3600;
TC_KF_config.init_clock_offset_unc = 10;
TC_KF_config.init_clock_drift_unc = 0.1;
TC_KF_config.init_clock_offset_unc_bds = 10;
TC_KF_config.init_clock_drift_unc_bds = 0.1;
             
TC_KF_config.gyro_noise_PSD = (0.01)^2;              
TC_KF_config.accel_noise_PSD = (0.1)^2;
TC_KF_config.accel_bias_PSD = 1.0E-5;
TC_KF_config.gyro_bias_PSD = 4.0E-11;
TC_KF_config.clock_freq_PSD = 1;
TC_KF_config.clock_phase_PSD = 1;  
TC_KF_config.pseudo_range_SD = 5;
TC_KF_config.range_rate_SD = 0.1;
omega_ie = 7.292115E-5;
Omega_ie = Skew_symmetric([0,0,omega_ie]); 
TC_KF_config.clock_freq_PSD_bds = 1;  
TC_KF_config.clock_phase_PSD_bds = 1;  
TC_KF_config.pseudo_range_SD_bds = 5;
TC_KF_config.range_rate_SD_bds = 0.1;

TC_KF_config.RecClockPreprocOptions = 2;

ARAIM_URA_GPS = 5; 
ARAIM_URA_GAL = 1;
ARAIM_URA_GLO = 1;
ARAIM_URA_BDU = 5; 

ARAIM_BIAS_GPS = .75;
ARAIM_BIAS_GAL = .75;
ARAIM_BIAS_GLO = .75;
ARAIM_BIAS_BDU = .75;

ARAIM_URE_GPS = 2/3*ARAIM_URA_GPS;
ARAIM_URE_GAL = 2/3*ARAIM_URA_GAL;
ARAIM_URE_GLO = 2/3*ARAIM_URA_GLO;
ARAIM_URE_BDU = 2/3*ARAIM_URA_BDU;

ARAIM_BIAS_CONT_GPS = 1;
ARAIM_BIAS_CONT_GAL = 1;
ARAIM_BIAS_CONT_GLO = 1;
ARAIM_BIAS_CONT_BDU = 1;

ARAIM_PSAT_GPS = 1e-5;
ARAIM_PSAT_GAL = 1e-5;
ARAIM_PSAT_GLO = 1e-3;
ARAIM_PSAT_BDU = 1e-5;

ARAIM_PCONST_GPS = 1e-8;
ARAIM_PCONST_GAL = 1e-4;
ARAIM_PCONST_GLO = 1e-4;
ARAIM_PCONST_BDU = 1e-8;

FDE_WF_FLAG = 0;

PL0_FDE = 0; 
ARAIM_USRMASK_GPS = mask_angle;
ARAIM_USRMASK_GLO = mask_angle;
ARAIM_USRMASK_GAL = mask_angle;
ARAIM_USRMASK_BDU = mask_angle;

ARAIM_SIN_USRMASK_GPS = sin(ARAIM_USRMASK_GPS*pi/180);
ARAIM_SIN_USRMASK_GLO = sin(ARAIM_USRMASK_GLO*pi/180);
ARAIM_SIN_USRMASK_GAL = sin(ARAIM_USRMASK_GAL*pi/180);
ARAIM_SIN_USRMASK_BDU = sin(ARAIM_USRMASK_BDU*pi/180);

PHMI_VERT = .1e-8; 
P_THRES = 6e-8;     
PFA_VERT = .01e-7;
PFA_HOR  = 1e-7;      
PL_TOL = 1e-2;
PHMI_HOR = 9.9e-8; 
P_EMT = 1e-5;
FC_THRES = .01;
 
VPLT = 35;
HPLT = 40;
EMTT = 15;

SIG_ACC_MAX_VERT = 10;
SIG_ACC_MAX_HOR1 = 10; 
SIG_ACC_MAX_HOR2 = 10;  

ATTEMPT_OPT = 0; 

disp(' ')
disp('...GNSS data loading...')
Rinex3_Parser
disp(' ')
disp('...IMU data loading...')
disp(' ')

first_epoch = 95593; 
last_epoch = 96379;  
original_obs_IMU = readmatrix("data/xsense_imu_medium_urban1.csv");
obs_IMU(:,1) = original_obs_IMU(:,1);
obs_IMU(:,2) = original_obs_IMU(:,3);
obs_IMU(:,3) = original_obs_IMU(:,2);
obs_IMU(:,4) = original_obs_IMU(:,19);
obs_IMU(:,5) = original_obs_IMU(:,18);
obs_IMU(:,6) = -original_obs_IMU(:,20);
obs_IMU(:,7) = original_obs_IMU(:,31);
obs_IMU(:,8) = original_obs_IMU(:,30);
obs_IMU(:,9) = -original_obs_IMU(:,32);
 

PRonL1_gps='C1C '; 
L1_gps='L1C ';     
D1_gps='D1C ';     
S1_gps='S1C ';     

PRonL1_bds='C2I ';
L1_bds='L2I ';    
D1_bds='D2I ';    
S1_bds='S2I ';    

Iono_corr_GPS = iono_corr_klo_gps;
Iono_corr_BDS = iono_corr_klo_bds;


observations.gps(:,2)=round(observations.gps(:,2));
observations.com(:,2)=round(observations.com(:,2));
epoche_obs=unique([observations.gps(:,2);observations.com(:,2)]);
nr_observations_gps=size(observations.gps,1);
nr_observations_bds=size(observations.com,1);
eph_GPS=eph.gps; 
eph_BDS=eph.com; 

epochs=first_epoch:dt_data:last_epoch;
epochs=makeitcol(epochs);                  

nr_epoche=length(epochs);
durata_sessione=last_epoch-first_epoch;

observations.gps=observations.gps(observations.gps(:,3)==0,:); 
obs_GPS=nan(nr_observations_gps,7); 
obs_GPS(:,1)=observations.gps(:,1); 
obs_GPS(:,2)=observations.gps(:,2);
obs_GPS(:,3)=observations.gps(:,4);
col_PRonL1_gps = find(ismember(header.gps,PRonL1_gps)==1);
if ~isempty(col_PRonL1_gps)
    obs_GPS(:,4)=observations.gps(:,col_PRonL1_gps);
end
col_L1_gps = find(ismember(header.gps,L1_gps)==1);
if ~isempty(col_L1_gps)
    obs_GPS(:,5)=observations.gps(:,col_L1_gps);
end
col_D1_gps = find(ismember(header.gps,D1_gps)==1);
if ~isempty(col_D1_gps)
    obs_GPS(:,6)=observations.gps(:,col_D1_gps);
end
col_S1_gps = find(ismember(header.gps,S1_gps)==1);
if ~isempty(col_S1_gps)
    obs_GPS(:,7)=observations.gps(:,col_S1_gps);
end

observations.com=observations.com(observations.com(:,3)==0,:);
obs_BDS=nan(nr_observations_bds,7);
obs_BDS(:,1)=observations.com(:,1);
obs_BDS(:,2)=observations.com(:,2);
obs_BDS(:,3)=observations.com(:,4);
col_C1_bds = find(ismember(header.com,PRonL1_bds)==1);
if ~isempty(col_C1_bds)
    obs_BDS(:,4)=observations.com(:,col_C1_bds);
end
col_L1_bds = find(ismember(header.com,L1_bds)==1);
if ~isempty(col_L1_bds)
    obs_BDS(:,5)=observations.com(:,col_L1_bds);
end
col_D1_bds = find(ismember(header.com,D1_bds)==1);
if ~isempty(col_D1_bds)
    obs_BDS(:,6)=observations.com(:,col_D1_bds);
end
col_S1_bds = find(ismember(header.com,S1_bds)==1);
if ~isempty(col_S1_bds)
    obs_BDS(:,7)=observations.com(:,col_S1_bds);
end

init_lat = deg2rad(22.301199858333334);
init_lon = deg2rad(114.1790571083333);
init_h = 3.472;
init_roll = deg2rad(0);
init_pitch = deg2rad(0);
init_yaw = deg2rad(-90);


[GNSS_init.x,GNSS_init.y,GNSS_init.z] = geo2ecef(init_lat,init_lon,init_h);
GNSS_init.vx = 0;    
GNSS_init.vy = 0;
GNSS_init.vz = 0;
GNSS_init.RecClk = 1e-4;
GNSS_init.RecClkDrift = 1e-6;
old_est_r_ea_e=[GNSS_init.x,GNSS_init.y,GNSS_init.z]'; 
old_est_v_ea_e=[GNSS_init.vx,GNSS_init.vy,GNSS_init.vz]';
est_clock(1:2)=[GNSS_init.RecClk*v_light, GNSS_init.RecClkDrift*v_light];
signalbds_init.RecClk = 0;
signalbds_init.RecClkDrift = 0;
est_clock(3:4)=[signalbds_init.RecClk*v_light, signalbds_init.RecClkDrift*v_light];


attitude_ini=[init_roll;init_pitch;init_yaw];
L_ba_b=[0;0.12;0.082]; 
old_est_C_b_n=Euler_to_CTM(attitude_ini)';
[old_est_L_a,old_est_lambda_a,old_est_h_a,old_est_v_ea_n] = pv_ECEF_to_NED(old_est_r_ea_e,old_est_v_ea_e); 
[~,~,old_est_C_b_e] = NED_to_ECEF(old_est_L_a,old_est_lambda_a,old_est_h_a,old_est_v_ea_n,old_est_C_b_n);   
old_est_r_eb_e=old_est_r_ea_e-old_est_C_b_e*L_ba_b; 
old_est_v_eb_e=old_est_v_ea_e; 
 
old_time = first_epoch;
Total_GNSS_epoch = nr_epoche;

init_const;      
init_col_labels; 
init_mops;       
[flag_profile,vhpl,gnssnum,out_profile,position_inner_std_profile,out_IMU_bias_est,out_clock,out_KF_SD,PickSubsetRes,RecClockBiasRes,...
    out_MeasurementNoise_SD,out_Resi,InnovationRes,StdInnovationRes] =...
    Tightly_coupled_INS_GNSS(eph_GPS,eph_BDS,obs_GPS,obs_BDS,obs_IMU,old_time,old_est_r_eb_e,old_est_v_eb_e,...
    est_clock,attitude_ini,GNSS_config,TC_KF_config,L_ba_b,Total_GNSS_epoch,...
    Iono_corr_GPS,Iono_corr_BDS,DOY);

