function [est_C_b_e_new,est_v_eb_e_new,est_r_eb_e_new,est_IMU_bias_new,...
            est_clock_new,P_matrix_new,R_matrix,H_matrix,K_matrix,W_matrix,C_matrix,V_matrix,Z_all,H_all] = TC_KF_Epoch(prn,...
            tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,...
            est_IMU_bias_old,est_clock_old,P_matrix_old,meas_f_ib_b,...
            TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
            G,delta_z,R,Z_all,H_all,GNSS_epoch)

global MOPS_MIN_GPSPRN MOPS_MAX_GPSPRN MOPS_MIN_GLOPRN MOPS_MAX_GLOPRN ...
       MOPS_MIN_GALPRN MOPS_MAX_GALPRN  ...
       MOPS_MIN_BDUPRN MOPS_MAX_BDUPRN 

global T_Global

omega_ie = 7.292115E-5;
R_0 = 6378137;
e = 0.0818191908425;
tau = 0.02;

prn_gps = find((prn >= MOPS_MIN_GPSPRN & prn <= MOPS_MAX_GPSPRN));
prn_gal = find((prn >= MOPS_MIN_GALPRN & prn <= MOPS_MAX_GALPRN));
prn_glo = find((prn >= MOPS_MIN_GLOPRN & prn <= MOPS_MAX_GLOPRN));
prn_bds = find((prn >= MOPS_MIN_BDUPRN & prn <= MOPS_MAX_BDUPRN));

nconst = sum([~isempty(prn_gps),~isempty(prn_bds)]);

[est_L_b_old,est_lambda_b_old,~,~,est_C_b_n_old] =...
    ECEF_to_NED(est_r_eb_e_old,est_v_eb_e_old,est_C_b_e_old);
% Skew symmetric matrix of Earth rate
Omega_ie = Skew_symmetric([0,0,omega_ie]);
       
Phi_matrix = eye(19);
Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) - Omega_ie * tor_s;
Phi_matrix(1:3,13:15) = est_C_b_e_old * tor_s;
Phi_matrix(4:6,1:3) = -tor_s * Skew_symmetric(est_C_b_e_old * meas_f_ib_b);
Phi_matrix(4:6,4:6) = Phi_matrix(4:6,4:6) - 2 * Omega_ie * tor_s;
geocentric_radius = R_0 / sqrt(1 - (e * sin(est_L_b_old))^2) *...
    sqrt(cos(est_L_b_old)^2 + (1 - e^2)^2 * sin(est_L_b_old)^2);
Phi_matrix(4:6,7:9) = -tor_s * 2 * Gravity_ECEF(est_r_eb_e_old) /...
    geocentric_radius * est_r_eb_e_old' / sqrt (est_r_eb_e_old' *...
    est_r_eb_e_old);
Phi_matrix(4:6,10:12) = est_C_b_e_old * tor_s;
Phi_matrix(7:9,4:6) = eye(3) * tor_s;
Phi_matrix(16,17) = tor_s;
Phi_matrix(18,19) = tor_s;

Q_prime_matrix = zeros(19);
Q_prime_matrix(1:3,1:3) = eye(3) * TC_KF_config.gyro_noise_PSD * tor_s;
Q_prime_matrix(4:6,4:6) = eye(3) * TC_KF_config.accel_noise_PSD * tor_s;
Q_prime_matrix(10:12,10:12) = eye(3) * TC_KF_config.accel_bias_PSD * tor_s;
Q_prime_matrix(13:15,13:15) = eye(3) * TC_KF_config.gyro_bias_PSD * tor_s;
Q_prime_matrix(16,16) = TC_KF_config.clock_phase_PSD * tor_s;
Q_prime_matrix(17,17) = TC_KF_config.clock_freq_PSD * tor_s;
Q_prime_matrix(18,18) = TC_KF_config.clock_phase_PSD_bds * tor_s;
Q_prime_matrix(19,19) = TC_KF_config.clock_freq_PSD_bds * tor_s;

x_est_propagated(1:19,1) = 0;
Clock_Reset_Flag = 1;
if Clock_Reset_Flag == 1 
    x_est_propagated(16,1)=0;
    x_est_propagated(17,1)=0;
    x_est_propagated(18,1)=0;
    x_est_propagated(19,1)=0;
else
    x_est_propagated(16,1) = est_clock_old(1) + est_clock_old(2) * tor_s;
    x_est_propagated(17,1) = est_clock_old(2);
    x_est_propagated(18,1) = est_clock_old(3) + est_clock_old(4) * tor_s;
    x_est_propagated(19,1) = est_clock_old(4);
end

P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) *...
    Phi_matrix' + 0.5 * Q_prime_matrix;

H_matrix = zeros(length(prn),19);
H_matrix(:,7:9) = G(:,1:3);
H_matrix(1:length(prn_gps),16) = ones(length(prn_gps),1);
H_matrix(length(prn_gps)+1:length(prn_gps)+length(prn_bds),18) = ones(length(prn_bds),1);

R_matrix = R;

for j = 1:length(prn_gps)
    delta_z(j) = delta_z(j) - x_est_propagated(16);       
end
for j = 1+length(prn_gps):length(prn_gps)+length(prn_bds) 
    delta_z(j) = delta_z(j) - x_est_propagated(18);
end

Z_all{GNSS_epoch,1} = delta_z;
H_all{GNSS_epoch} = H_matrix;

if isempty(Z_all{GNSS_epoch-1}) || any(size(Z_all{GNSS_epoch-1}) == 0)
    historical_residual = delta_z;
else
    historical_residual = Z_all{GNSS_epoch-1};
end

len_delta_z = length(delta_z);
len_historical_residual = length(historical_residual);

if len_delta_z > len_historical_residual
    if len_historical_residual >= 2
        historical_residual_aligned = interp1(1:len_historical_residual, historical_residual, ...
                                              linspace(1, len_historical_residual, len_delta_z));
    else
        historical_residual_aligned = repmat(historical_residual, len_delta_z, 1);
    end
    delta_z_aligned = delta_z;
else
    if len_delta_z >= 2
        delta_z_aligned = interp1(1:len_delta_z, delta_z, ...
                                  linspace(1, len_delta_z, len_historical_residual));
    else
        delta_z_aligned = repmat(delta_z, len_historical_residual, 1);
    end
    historical_residual_aligned = historical_residual;
end

if norm(historical_residual_aligned) < eps
    residual_change_rate = 1e-9;
else
    residual_change_rate = norm(delta_z_aligned - historical_residual_aligned) / norm(historical_residual_aligned);
end

phi_min = 0.001;
phi_max = 1;

phi_z = tau * (1 + residual_change_rate);
phi_z = max(phi_min, min(phi_z, phi_max));

size1 = length(prn);

if isempty(Z_all{GNSS_epoch-1}) || any(size(Z_all{GNSS_epoch-1}) == 0)
    Z_k_star = delta_z;
    H_k_star = H_matrix;
else
    phi_v = ones(size1, size(Z_all{GNSS_epoch-1}, 1));
    phi = phi_z/size(Z_all{GNSS_epoch-1}, 1);
    Z_k_star = delta_z - phi * phi_v * Z_all{GNSS_epoch-1, 1};
    H_k_star = H_matrix * Phi_matrix - phi * phi_v * H_all{GNSS_epoch-1};
end

R_k_star = H_matrix * Q_prime_matrix * H_matrix' + R_matrix;
W = H_k_star * P_matrix_propagated * H_k_star' + R_k_star;

Y_matrix = [Z_k_star;x_est_propagated];
C_matrix = [H_k_star;eye(19)];
[size1,size2] = size(R_matrix);
V_matrix = [R_k_star,zeros(size1,19);zeros(19,size2),P_matrix_propagated];
W_matrix = pinv(V_matrix);
x_1 = C_matrix' * W_matrix * C_matrix;
x_2 = pinv(x_1);
x_est_new = x_2 * C_matrix' * W_matrix * Y_matrix;

K_matrix = P_matrix_propagated * H_k_star' / W;

P_matrix_new = (eye(19) - K_matrix * H_k_star) * P_matrix_propagated;

est_C_b_e_new = (eye(3) - Skew_symmetric(x_est_new(1:3))) * est_C_b_e_old;
est_v_eb_e_new = est_v_eb_e_old - x_est_new(4:6);
est_r_eb_e_new = est_r_eb_e_old - x_est_new(7:9);

est_IMU_bias_new = est_IMU_bias_old + x_est_new(10:15);
est_clock_new = x_est_new(16:19)';

% Ends

