function [sigma, sigma_ss, bias, bias_ss, s1vec, s2vec, s3vec, x, So,...
    lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,Z_all,H_all] = compute_subset_EKF2_solutions(...
             G, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, subsets, delta_z,...
             tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,est_IMU_bias_old,est_clock_old,P_matrix_old,...
                            meas_f_ib_b_old,TC_KF_config,L_ba_b,meas_omega_ib_b_old,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,R,...
                            prn, el, snr,Z_all,H_all,GNSS_epoch)

N=size(G,1);
nsets=size(subsets,1);

lamda = 0;

s1vec = ones(nsets,N)*Inf;
s2vec = ones(nsets,N)*Inf;
s3vec = ones(nsets,N)*Inf;
dP1 = ones(nsets,1)*Inf;
dP2 = ones(nsets,1)*Inf;
dP3 = ones(nsets,1)*Inf;
x     = ones(nsets,size(G,2));

for i=1:nsets      
    index = find(subsets(i,:)>0);
    index0 = find(subsets(i,:)==0);

    if (~isempty(index))  
        
        [est_C_b_e_s,est_v_eb_e_s,est_r_eb_e_s,est_IMU_bias_s,est_clock_s,P_matrix_s,~,H_matrix_s,K_matrix_s,W_matrix_s,C_matrix_s,V_matrix_s,Z_all,H_all] =...
                my_TC_KF_Epoch_8_2(prn(index),tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,est_IMU_bias_old,est_clock_old,P_matrix_old,...
                                meas_f_ib_b_old,TC_KF_config,L_ba_b,meas_omega_ib_b_old,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
                                G(index,:),delta_z(index,:),R(index,index),Z_all,H_all,GNSS_epoch);
  
        if (i==1)
            C = C_matrix_s;
            W = W_matrix_s;
            CtW = C'*W;
            invcov0 = CtW*C;
            cov0 = inv(invcov0);
        end
        S = compute_s_coefficients(C, W, CtW, invcov0, index0, cov0);

        if(i==1)
            So = S; 
            est_C_b_e = est_C_b_e_s;
            est_v_eb_e = est_v_eb_e_s;
            est_r_eb_e = est_r_eb_e_s;
            est_IMU_bias = est_IMU_bias_s;
            est_clock = est_clock_s;
            P_matrix = P_matrix_s;
        end
        P = P_matrix_s(7:9,7:9);
        dP = P - P_matrix(7:9,7:9);

    else
        So = 0;
        S = ones(3,N)*Inf;
        est_r_eb_e_s = ones(3,1)*Inf;
        est_clock_s = ones(1,6)*Inf;
    end   

        s1vec(i,:)= S(7,1:length(sigpr2_int));
        s2vec(i,:)= S(8,1:length(sigpr2_int));
        s3vec(i,:)= S(9,1:length(sigpr2_int));
        dP1(i) = sqrt(abs(dP(1,1)));
        dP2(i) = sqrt(abs(dP(2,2)));
        dP3(i) = sqrt(abs(dP(3,3)));
        x(i,1:3) = est_r_eb_e_s';
        x(i,4) = est_clock_s(1);
        x(i,5) = est_clock_s(3);
end

V = V_matrix_s(1:length(sigpr2_int),1:length(sigpr2_int));

sigma(:,1) = dP1;
sigma(:,2) = dP2;
sigma(:,3) = dP3;

bias(:,1)   = abs(s1vec)*nom_bias_int;
bias(:,2)   = abs(s2vec)*nom_bias_int;
bias(:,3)   = abs(s3vec)*nom_bias_int;

delta_s1vec = s1vec - ones(nsets,1)*s1vec(1,:);
delta_s2vec = s2vec - ones(nsets,1)*s2vec(1,:);
delta_s3vec = s3vec - ones(nsets,1)*s3vec(1,:);

sigma_ss(:,1) = sqrt(diag(delta_s1vec*V*delta_s1vec'));
sigma_ss(:,2) = sqrt(diag(delta_s2vec*V*delta_s2vec'));
sigma_ss(:,3) = sqrt(diag(delta_s3vec*V*delta_s3vec'));

bias_ss(:,1)   = abs(delta_s1vec)*nom_bias_acc;
bias_ss(:,2)   = abs(delta_s2vec)*nom_bias_acc;
bias_ss(:,3)   = abs(delta_s3vec)*nom_bias_acc;

%End 



