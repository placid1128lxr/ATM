function [FLAG_Fault,return_lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
    vpl, hpl, sig_acc, emt, subsets, pap_subset, p_not_monitored, rho,Z_all,H_all] = mhss_raim_baseline(prn,tor_s,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                            meas_f_ib_b,TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,R,...
                            delta_zz,G, el, snr,sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,...
    opt_flag, rho_j,Z_all,H_all,GNSS_epoch, subsets, pap_subset, p_not_monitored)

rho = [];

if nargin<27
    rho_j = 1;
end

global PHMI_VERT PHMI_HOR P_THRES PFA_VERT PFA_HOR P_EMT PL_TOL FC_THRES
global SIG_ACC_MAX_VERT SIG_ACC_MAX_HOR1 SIG_ACC_MAX_HOR2 

%%%%%%%%%%%% 1 Determine subsets and associated probabilities %%%%%%%%%%%%%%%
if nargin<30 
 [subsets, pap_subset, p_not_monitored] = determine_subsets_v4(G, p_sat, p_const, P_THRES, FC_THRES);
end


%%%%%%%%%%%%% 2 Compute subset position solutions, sigmas, and biases %%%%%%%
[sigma, sigma_ss, bias, bias_ss, s1vec, s2vec, s3vec, x, ~,...
    ~,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,Z_all,H_all] = compute_subset_EKF2_solutions(...
             G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, subsets, delta_zz,...
             tor_s,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                            meas_f_ib_b,TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,R,...
                            prn, el, snr,Z_all,H_all,GNSS_epoch);

%%%%%%%%%%%%% 3 Adjust all-in-view position coefficients %%%%%%%%%%%%%%%%%%%%
if opt_flag
nsets=size(subsets,1);

[ s3opt ] = compute_adjusted_position_1D(pap_subset,sigma(:,3),sigma_ss(:,3), s3vec, sigpr2_acc,SIG_ACC_MAX_VERT);
[ s2opt ] = compute_adjusted_position_1D(pap_subset,sigma(:,2),sigma_ss(:,2), s2vec, sigpr2_acc,SIG_ACC_MAX_HOR2);
[ s1opt ] = compute_adjusted_position_1D(pap_subset,sigma(:,1),sigma_ss(:,1), s1vec, sigpr2_acc,SIG_ACC_MAX_HOR1);

s1vec(1,:) = s1opt;  
s2vec(1,:) = s2opt;  
s3vec(1,:) = s3opt;  
sigma(1,1) = sqrt((s1vec(1,:).^2)* sigpr2_int);
sigma(1,2) = sqrt((s2vec(1,:).^2)* sigpr2_int);
sigma(1,3) = sqrt((s3vec(1,:).^2)* sigpr2_int);

bias(1,1)   = abs(s1vec(1,:))*nom_bias_int;
bias(1,2)   = abs(s2vec(1,:))*nom_bias_int;
bias(1,3)   = abs(s3vec(1,:))*nom_bias_int;

delta_s1vec = s1vec - ones(nsets,1)*s1vec(1,:);
delta_s2vec = s2vec - ones(nsets,1)*s2vec(1,:);
delta_s3vec = s3vec - ones(nsets,1)*s3vec(1,:);

sigma_ss(:,1) = sqrt((delta_s1vec.^2)* sigpr2_acc);
sigma_ss(:,2) = sqrt((delta_s2vec.^2)* sigpr2_acc);
sigma_ss(:,3) = sqrt((delta_s3vec.^2)* sigpr2_acc);

bias_ss(:,1)   = abs(delta_s1vec)*nom_bias_acc;
bias_ss(:,2)   = abs(delta_s2vec)*nom_bias_acc;
bias_ss(:,3)   = abs(delta_s3vec)*nom_bias_acc;

end        


%%%%%%%%%%%%% 4 Filter out modes that cannot be monitored and adjust%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% integrity budget %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [sigma, sigma_ss, bias, bias_ss, s1vec, s2vec, s3vec, subsets,...
     pap_subset, p_not_monitored_2, ~, x] = filter_out_subsets(sigma, sigma_ss, bias,...
     bias_ss, s1vec, s2vec, s3vec, subsets, pap_subset, p_not_monitored, x);
 
 if nargin<30
     added_p_not_mon = 0;
     p_not_monitored = p_not_monitored_2;
 else
     added_p_not_mon = p_not_monitored_2 - p_not_monitored;
 end

%%%%%%%%%%%%%%%%%%%%%%% 5 Compute test thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%
[T1, T2, T3] = compute_test_thresholds(sigma_ss, bias_ss, PFA_VERT, PFA_HOR);

nsets=size(pap_subset,1);
delta_xx = x(:,1) - ones(nsets,1)*x(1,1);
delta_yy = x(:,2) - ones(nsets,1)*x(1,2);
delta_zz = x(:,3) - ones(nsets,1)*x(1,3);
FLAG_Fault = 0;
if ( find((abs(delta_xx)-T1)>0) )
    %disp('Have Fault!')
    FLAG_Fault = 1;
end
if ( find((abs(delta_yy)-T2)>0) )
    %disp('Have Fault!')
    FLAG_Fault = 1;
end
if ( find((abs(delta_zz)-T3)>0) )
    %disp('Have Fault!')
    FLAG_Fault = 1;
end

lamda1 = sqrt(delta_xx.^2 + delta_yy.^2)./sqrt(sigma_ss(:,1)+sigma_ss(:,2));
lamda2 = sqrt(delta_xx.^2 + delta_yy.^2 + delta_zz.^2)./sqrt(sigma_ss(:,1)+sigma_ss(:,2)+sigma_ss(:,3));
lamda3 = sqrt(delta_xx.^2 + delta_yy.^2);
lamda4 = sqrt(delta_xx.^2 + delta_yy.^2 + delta_zz.^2);
lamda5 = sqrt(delta_xx.^2)./sqrt(sigma_ss(:,1)) + sqrt(delta_yy.^2)./sqrt(sigma_ss(:,2)) + sqrt(delta_zz.^2)./sqrt(sigma_ss(:,3));
lamda6 = sqrt(1*delta_xx.^2 + 1*delta_yy.^2 + 0.1*delta_zz.^2);
lamda = lamda4;

return_lamda = ones(length(prn),1);
for i = 2:nsets 
    lamda_idx = find(subsets(i,:)==0);
    if( length(lamda_idx)==1 )
        return_lamda(lamda_idx) = lamda(i);
    end
end

%%%%%%%%%%%%%%%%%%%%% 6 Compute Vertical Protection Level %%%%%%%%%%%%%%%%%%%

p_fault = pap_subset;
p_fault(1) = 2;
phmi_vert = rho_j*(PHMI_VERT - PHMI_VERT/(PHMI_VERT+PHMI_HOR)*p_not_monitored) -PHMI_VERT/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[vpl, alloc3]     = compute_protection_level_v4_2(sigma(:,3), bias(:,3) + T3, p_fault, phmi_vert, PL_TOL);


%%%%%%%%%%%%%%%%%%%%% 6 Compute Horizontal Protection Level %%%%%%%%%%%%%%%%%
phmi_hor = rho_j*(PHMI_HOR/2 -.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*p_not_monitored)-.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[hpl1, alloc1] = compute_protection_level_v4_2(sigma(:,1), bias(:,1) + T1, p_fault,  phmi_hor, PL_TOL);
[hpl2, alloc2] = compute_protection_level_v4_2(sigma(:,2), bias(:,2) + T2, p_fault,  phmi_hor, PL_TOL);
hpl = sqrt(hpl1^2 + hpl2^2);

%7 Exclude modes that are double counted

idx = find((alloc3+alloc2+alloc1)>=1);

if ~isempty(idx)

p_not_monitored = p_not_monitored + sum(p_fault(idx));
idx = setdiff(1:length(p_fault),idx);
s1vec = s1vec(idx,:);
s2vec = s2vec(idx,:);
s3vec = s3vec(idx,:);
sigma = sigma(idx,:);
sigma_ss = sigma_ss(idx,:);      
bias = bias(idx,:);
bias_ss = bias_ss(idx,:);
subsets = subsets(idx,:);
%p_fault =p_fault(idx);
pap_subset= pap_subset(idx);
%Nsubsets =length(idx);

%%%%%%%%%%%%%%%%%%%%%%% 8 Compute test thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%

[T1, T2, T3] = compute_test_thresholds(sigma_ss, bias_ss, PFA_VERT, PFA_HOR);

%%%%%%%%%%%%%%%%%%%%% 9 Compute Vertical Protection Level %%%%%%%%%%%%%%%%%%%

p_fault = pap_subset;
p_fault(1) = 2;
phmi_vert = rho_j*(PHMI_VERT - PHMI_VERT/(PHMI_VERT+PHMI_HOR)*p_not_monitored) -PHMI_VERT/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[vpl, alloc3]     = compute_protection_level_v4_2(sigma(:,3), bias(:,3) + T3, p_fault,  phmi_vert, PL_TOL);

%%%%%%%%%%%%%%%%%%%%% 9 Compute Horizontal Protection Level %%%%%%%%%%%%%%%%%
phmi_hor = rho_j*(PHMI_HOR/2 -.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*p_not_monitored)-.5*PHMI_HOR/(PHMI_VERT+PHMI_HOR)*added_p_not_mon ;
[hpl1, alloc1] = compute_protection_level_v4_2(sigma(:,1), bias(:,1) + T1, p_fault,  phmi_hor, PL_TOL,1-alloc3);
[hpl2, ~] = compute_protection_level_v4_2(sigma(:,2), bias(:,2) + T2, p_fault,  phmi_hor, PL_TOL, 1-alloc1-alloc3);
hpl = sqrt(hpl1^2 + hpl2^2);


end


%%%%%%%%%%%%%%%%%%% 10 Compute Effective Monitor threshold %%%%%%%%%%%%%%%%%%%
% idx = find(pap_subset>=P_EMT);
sigma3_acc = sqrt((s3vec.*s3vec)*sigpr2_acc);
% K_md_emt   = - norminv(.5*P_EMT./pap_subset);
%emt = max(T3(idx) + K_md_emt(idx).*sigma3_acc(idx));
emt = compute_emt(s3vec, sigpr2_acc, T3, pap_subset, P_EMT);

%%%%%%%%%%%%%%%%%%%% 11 Compute sigma accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_acc = sigma3_acc(1);
%End





