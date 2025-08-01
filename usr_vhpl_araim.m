function [vhpl,lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,Z_all,H_all]=usr_vhpl_araim(usr_idx, sig2_i, prn, sig2acc_i, bnom_i, bcont_i, psat_i,...
    GNSS_measurements,...
            no_meas,tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,...
            est_IMU_bias_old,est_clock_old,P_matrix_old,meas_f_ib_b,...
            TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
            sigma_2old,Z_all,H_all,GNSS_epoch)

global MOPS_NOT_MONITORED 

global MOPS_MIN_GPSPRN MOPS_MAX_GPSPRN MOPS_MIN_GLOPRN MOPS_MAX_GLOPRN ...
       MOPS_MIN_GALPRN MOPS_MAX_GALPRN  ...
       MOPS_MIN_BDUPRN MOPS_MAX_BDUPRN                                         
global ARAIM_PCONST_GPS ARAIM_PCONST_GAL ARAIM_PCONST_GLO ARAIM_PCONST_BDU  
global SIG_ACC_MAX_VERT SIG_ACC_MAX_HOR1 SIG_ACC_MAX_HOR2 VPLT HPLT EMTT FDE_FLAG
global ATTEMPT_OPT PL0_FDE FDE_WF_FLAG

global init_lat init_lon init_h  
global flag_fault_record  flag_fault_record_bl

global PHMI_VERT PHMI_HOR P_THRES PFA_VERT PFA_HOR P_EMT PL_TOL FC_THRES

global flag_opt T_Global

global mhss_fault GDOP

mhss_fault = 0;
GDOP = Inf;

%initialize return value
n_usr=max(usr_idx);

vhpl=repmat(MOPS_NOT_MONITORED,n_usr,4);

for usr=1:n_usr
    
    svidxgps = find((usr_idx==usr)&(prn >= MOPS_MIN_GPSPRN & prn <= MOPS_MAX_GPSPRN)&(sig2_i<Inf)&(psat_i<1));
    svidxgal = find((usr_idx==usr)&(prn >= MOPS_MIN_GALPRN & prn <= MOPS_MAX_GALPRN)&(sig2_i<Inf)&(psat_i<1));
    svidxglo = find((usr_idx==usr)&(prn >= MOPS_MIN_GLOPRN & prn <= MOPS_MAX_GLOPRN)&(sig2_i<Inf)&(psat_i<1));
    svidxbdu = find((usr_idx==usr)&(prn >= MOPS_MIN_BDUPRN & prn <= MOPS_MAX_BDUPRN)&(sig2_i<Inf)&(psat_i<1));

    ngps = length(svidxgps);
    ngal = length(svidxgal);
    nglo = length(svidxglo);
    nbdu = length(svidxbdu);
  
    svconst = [(~isempty(svidxgps)) (~isempty(svidxgal)) (~isempty(svidxglo)) (~isempty(svidxbdu))]; 
  
    nsat   = ngps + ngal + nglo + nbdu; 
    svidx  = [svidxgps; svidxgal; svidxglo; svidxbdu];

    clk_gps = [ones(ngps,svconst(1)); zeros(nsat - ngps,svconst(1))];
    clk_gal = [zeros(ngps,svconst(2)); ones(ngal,svconst(2));zeros(nsat-ngps-ngal,svconst(2))];    
    clk_glo = [zeros(ngps+ngal,svconst(3)); ones(nglo,svconst(3));zeros(nsat-ngps-ngal-nglo,svconst(3))];
    clk_bdu = [zeros(ngps+ngal+nglo,svconst(4)); ones(nbdu,svconst(4))];

    est_r_ea_e_old = est_r_eb_e_old + est_C_b_e_old*L_ba_b;
    u_as_e_T = zeros(no_meas,3);
    pred_meas = zeros(no_meas,1);    

    for j = 1:no_meas
        delta_r =  GNSS_measurements(j,3:5)' - est_r_ea_e_old;
        range = sqrt(delta_r' * delta_r); 
        pred_meas(j,1) = range;
        u_as_e_T(j,1:3) = delta_r' / range; 
    end 

 
    los_xyzb = zeros(size(GNSS_measurements,1),3);
    los_xyzb(:,1:3) = u_as_e_T(1:no_meas,1:3);
    G = [ los_xyzb(svidx,1:3) clk_gps clk_gal clk_glo clk_bdu]; 

 
    delta_z = zeros(size(GNSS_measurements,1),1);
    delta_z(1:no_meas,1) = GNSS_measurements(1:no_meas,1) - pred_meas(1:no_meas,1);
    el = GNSS_measurements(:,9);
    snr = GNSS_measurements(:,10);
 
    ElevationGain_GPS = (1./(sind(el(1:ngps))).^2);
    ElevationGain_BDS = (1./(sind(el(1+ngps:ngps+nbdu))).^2);
  
    R_matrix(1:ngps,1:ngps) = diag(ElevationGain_GPS)*...
    TC_KF_config.pseudo_range_SD^2;
    R_matrix(1+ngps:ngps+nbdu,1+ngps:ngps+nbdu) = ...
         diag(ElevationGain_BDS)*TC_KF_config.pseudo_range_SD_bds^2;


    p_const = [ones(svconst(1),1)*ARAIM_PCONST_GPS; ones(svconst(2),1)*ARAIM_PCONST_GAL;...
               ones(svconst(3),1)*ARAIM_PCONST_GLO; ones(svconst(4),1)*ARAIM_PCONST_BDU];
   
    sigpr2_int = sig2_i(svidx);
    sigpr2_acc = sig2acc_i(svidx);
    nom_bias_int = bnom_i(svidx);
    nom_bias_acc = bcont_i(svidx);
    p_sat = psat_i(svidx);

    n_view = length(svidx);
  
    if(n_view > 3 && FDE_FLAG == 1 )
        
            est_C_b_e_for_EKF = est_C_b_e_old;
            est_v_eb_e_for_EKF = est_v_eb_e_old;
            est_r_eb_e_for_EKF = est_r_eb_e_old;
            est_IMU_bias_for_EKF = est_IMU_bias_old;
            est_clock_for_EKF = est_clock_old;
            P_matrix_for_EKF = P_matrix_old;

                opt_idx = true(1,n_view);

                svidx = svidx(opt_idx);
                n_view = length(svidx);

            nn_view = n_view;
            if(nn_view > 3) 

                prn = prn(opt_idx);
                R_matrix = R_matrix(opt_idx,opt_idx);
                delta_z = delta_z(opt_idx);
                el = el(opt_idx);
                snr = snr(opt_idx);
                sigpr2_int = sigpr2_int(opt_idx); 
                sigpr2_acc = sigpr2_acc(opt_idx);
                nom_bias_int = nom_bias_int(opt_idx);
                nom_bias_acc =nom_bias_acc(opt_idx); 
                p_sat = p_sat(opt_idx);

                idx_G = find(sum(G(opt_idx,4:end))==0)+3;
                opt_G = G(opt_idx,:);
                opt_G(:,idx_G) = [];
                G = opt_G;
                subsets_exc = ones(nn_view,nn_view) - eye(nn_view);

                Nsubsets = size(subsets_exc,1);
                Nexc = size(subsets_exc,1);
                hpl_exc = Inf*ones(Nsubsets,1);
                vpl_exc = Inf*ones(Nsubsets,1);
                sig_acc_exc = Inf*ones(Nsubsets,1);
                emt_exc     = Inf*ones(Nsubsets,1);
                rho_j = 1/(Nexc+1);

                [subsets, pap_subset, p_not_monitored] = determine_subsets_v4(G, p_sat, p_const, P_THRES, FC_THRES);

                [FLAG_Fault,lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                    vpl, hpl, sig_acc, emt,subsets, pap_subset, p_not_monitored,~,Z_all,H_all] = ...
                    mhss_raim_baseline(prn,...
                    tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,...
                    est_IMU_bias_old,est_clock_old,P_matrix_old,meas_f_ib_b,...
                    TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
                    R_matrix,...
                    delta_z, G, el, snr, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,0, rho_j,Z_all,H_all,GNSS_epoch,subsets, pap_subset, p_not_monitored);

                for i=1:Nsubsets
                    idx = logical(subsets_exc(i,:));
                    subsets_i = subsets(:,idx);
                    if sum(idx)>3

                        [FLAG_Fault_bl,lamda_bl,est_C_b_e_bl,est_v_eb_e_bl,est_r_eb_e_bl,est_IMU_bias_bl,est_clock_bl,P_matrix_bl,...
                            vpl_bl, hpl_bl, sig_acc_bl, emt_bl,~,~,~,~,Z_all,H_all] = ...
                            mhss_raim_baseline(prn(idx),...
                            tor_s,est_C_b_e_for_EKF,est_v_eb_e_for_EKF,est_r_eb_e_for_EKF,...
                            est_IMU_bias_for_EKF,est_clock_for_EKF,P_matrix_for_EKF,meas_f_ib_b,...
                            TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
                            R_matrix(idx,idx),...
                            delta_z(idx), G(idx,:), el(idx), snr(idx), ...
                            sigpr2_int(idx), sigpr2_acc(idx),nom_bias_int(idx), nom_bias_acc(idx), p_sat(idx), p_const, 0, rho_j, Z_all,H_all,GNSS_epoch,...
                            subsets_i, pap_subset, p_not_monitored);

                        vpl_exc(i) = vpl_bl;
                        emt_exc(i) = emt_bl;
                        sig_acc_exc(i) = sig_acc_bl;
                        hpl_exc(i) = hpl_bl;

                        FLAG_Fault_bl_record(i) = FLAG_Fault_bl;
                        est_C_b_e_bl_record{i} = est_C_b_e_bl;
                        est_v_eb_e_bl_record{i} = est_v_eb_e_bl;
                        est_r_eb_e_bl_record{i} = est_r_eb_e_bl;
                        est_IMU_bias_bl_record{i} = est_IMU_bias_bl;
                        est_clock_bl_record{i} = est_clock_bl;
                        P_matrix_bl_record{i} = P_matrix_bl;
                        lamda_bl_record{i} = lamda_bl;

                    else
                        FLAG_Fault_bl_record = [];
                    end

                end


                hpl = max(hpl,max(hpl_exc(1:Nsubsets)));
                vpl = max(vpl,max(vpl_exc(1:Nsubsets)));
                emt = max(emt,max(emt_exc(1:Nsubsets)));
                sig_acc = max(sig_acc, max(sig_acc_exc(1:Nsubsets)));

                if(FLAG_Fault == 1)
                    idx_without_fault = find(FLAG_Fault_bl_record==0);
                    if (~isempty(idx_without_fault))
                        est_C_b_e = est_C_b_e_bl_record{idx_without_fault(1)};
                        est_v_eb_e = est_v_eb_e_bl_record{idx_without_fault(1)};
                        est_r_eb_e = est_r_eb_e_bl_record{idx_without_fault(1)};
                        est_IMU_bias = est_IMU_bias_bl_record{idx_without_fault(1)};
                        est_clock = est_clock_bl_record{idx_without_fault(1)};
                        P_matrix = P_matrix_bl_record{idx_without_fault(1)};
                        mhss_fault = 1;


                        emt = emt_exc(idx_without_fault(1));
                        vpl = vpl_exc(idx_without_fault(1));
                        sig_acc = sig_acc_exc(idx_without_fault(1));
                        hpl = hpl_exc(idx_without_fault(1));

                        sub_idx = logical(subsets_exc(idx_without_fault(1),:));
                        svidx = svidx(sub_idx);
                        n_view = length(svidx);

                        nn_view = n_view;
                        subsets_exc = ones(nn_view,nn_view) - eye(nn_view);

                        Nsubsets = size(subsets_exc,1);
                        Nexc = size(subsets_exc,1);
                        hpl_exc = Inf*ones(Nsubsets,1);
                        vpl_exc = Inf*ones(Nsubsets,1);
                        sig_acc_exc = Inf*ones(Nsubsets,1);
                        emt_exc     = Inf*ones(Nsubsets,1);
                        rho_j = 1/(Nexc+1);

                        for i=1:Nsubsets
                            idx = logical(subsets_exc(i,:));
                            subsets_i = subsets(:,idx);
                            if sum(idx)>3

                                [FLAG_Fault_bl,lamda_bl,est_C_b_e_bl,est_v_eb_e_bl,est_r_eb_e_bl,est_IMU_bias_bl,est_clock_bl,P_matrix_bl,...
                                    vpl_bl, hpl_bl, sig_acc_bl, emt_bl,~,~,~,~,Z_all,H_all] = ...
                                    mhss_raim_baseline(prn(idx),...
                                    tor_s,est_C_b_e_for_EKF,est_v_eb_e_for_EKF,est_r_eb_e_for_EKF,...
                                    est_IMU_bias_for_EKF,est_clock_for_EKF,P_matrix_for_EKF,meas_f_ib_b,...
                                    TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
                                    R_matrix(idx,idx),...
                                    delta_z(idx), G(idx,:), el(idx), snr(idx), ...
                                    sigpr2_int(idx), sigpr2_acc(idx),nom_bias_int(idx), nom_bias_acc(idx), p_sat(idx), p_const, 0, rho_j,Z_all,H_all,GNSS_epoch, ...
                                    subsets_i, pap_subset, p_not_monitored);

                                vpl_exc(i) = vpl_bl;
                                emt_exc(i) = emt_bl;
                                sig_acc_exc(i) = sig_acc_bl;
                                hpl_exc(i) = hpl_bl;

                                FLAG_Fault_bl_record(i) = FLAG_Fault_bl;
                                est_C_b_e_bl_record{i} = est_C_b_e_bl;
                                est_v_eb_e_bl_record{i} = est_v_eb_e_bl;
                                est_r_eb_e_bl_record{i} = est_r_eb_e_bl;
                                est_IMU_bias_bl_record{i} = est_IMU_bias_bl;
                                est_clock_bl_record{i} = est_clock_bl;
                                P_matrix_bl_record{i} = P_matrix_bl;

                            else
                                FLAG_Fault_bl_record = [];
                            end

                        end

                        hpl = max(hpl,max(hpl_exc(1:Nsubsets)));
                        vpl = max(vpl,max(vpl_exc(1:Nsubsets)));
                        emt = max(emt,max(emt_exc(1:Nsubsets)));
                        sig_acc = max(sig_acc, max(sig_acc_exc(1:Nsubsets)));
                        
                    else 
                        mhss_fault = 2;
                        disp('Still Have Fault!')
                    end

                end

            else

                if PL0_FDE
                    rho_j = 1/n_view;
                else
                    rho_j = 1;
                end
                [FLAG_Fault,lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,vpl, hpl, sig_acc, emt,~,~,~,~,Z_all,H_all] = ...
                    mhss_raim_baseline(prn,...
                    tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,...
                    est_IMU_bias_old,est_clock_old,P_matrix_old,meas_f_ib_b,...
                    TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
                    R_matrix,...
                    delta_z, G, el, snr, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,0, rho_j,Z_all,H_all,GNSS_epoch,subsets, pap_subset, p_not_monitored);
            end

            vhpl(usr,1)= vpl;
            vhpl(usr,2)= emt;
            vhpl(usr,3)= sig_acc;
            vhpl(usr,4)= hpl;
    else      
        if PL0_FDE
            rho_j = 1/n_view;
        else
            rho_j = 1;
        end 

        [subsets, pap_subset, p_not_monitored] = determine_subsets_v4(G, p_sat, p_const, P_THRES, FC_THRES);
        [FLAG_Fault,lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,vpl, hpl, sig_acc, emt,~,~,~,~,Z_all,H_all] = ...
            mhss_raim_baseline(prn,...
            tor_s,est_C_b_e_old,est_v_eb_e_old,est_r_eb_e_old,...
            est_IMU_bias_old,est_clock_old,P_matrix_old,meas_f_ib_b,...
            TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,...
            R_matrix,...
            delta_z, G, el, snr, sigpr2_int, sigpr2_acc, nom_bias_int, nom_bias_acc, p_sat, p_const,0, rho_j,Z_all,H_all,GNSS_epoch,subsets, pap_subset, p_not_monitored);


       
        vhpl(usr,1)=Inf;
        vhpl(usr,2)=Inf;
        vhpl(usr,3)=Inf;
        vhpl(usr,4)=Inf;
    
    end
end

