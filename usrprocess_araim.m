function [vhpl,lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,Z_all,H_all] = usrprocess_araim(prn_gnss,...
            el_s,az_s,...
            measures_all,no_GNSS_meas,tor_s,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                            meas_f_ib_b,TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,lamda,Z_all,H_all,GNSS_epoch)

global ARAIM_URA_GPS ARAIM_URA_GLO ARAIM_URA_GAL ARAIM_URA_BDU ...
       ARAIM_BIAS_GPS ARAIM_BIAS_GLO ARAIM_BIAS_GAL ARAIM_BIAS_BDU

global ARAIM_URE_GPS ARAIM_URE_GLO ARAIM_URE_GAL ARAIM_URE_BDU ...
       ARAIM_BIAS_CONT_GPS ARAIM_BIAS_CONT_GLO ARAIM_BIAS_CONT_GAL ARAIM_BIAS_CONT_BDU

global ARAIM_PSAT_GPS ARAIM_PSAT_GAL ARAIM_PSAT_GLO ARAIM_PSAT_BDU 

global MOPS_MIN_GPSPRN MOPS_MAX_GPSPRN MOPS_MIN_GLOPRN MOPS_MAX_GLOPRN ...
       MOPS_MIN_GALPRN MOPS_MAX_GALPRN MOPS_MIN_BDUPRN MOPS_MAX_BDUPRN

global ARAIM_USRMASK_GPS ARAIM_USRMASK_GLO ARAIM_USRMASK_GAL ARAIM_USRMASK_BDU 

global GAMMA_L1_L5 CONST_GAMMA_L1L5

global flag_sig2_if

global prn_index

    prn = prn_gnss;
    idusrgps = find(prn >= MOPS_MIN_GPSPRN & prn <= MOPS_MAX_GPSPRN);
    idusrglo = find(prn >= MOPS_MIN_GLOPRN & prn <= MOPS_MAX_GLOPRN);
    idusrgal = find(prn >= MOPS_MIN_GALPRN & prn <= MOPS_MAX_GALPRN);
    idusrbdu = find(prn >= MOPS_MIN_BDUPRN & prn <= MOPS_MAX_BDUPRN);
    
    abv_mask_gps = measures_all(idusrgps,9) >= ARAIM_USRMASK_GPS;
    abv_mask_glo = measures_all(idusrglo,9) >= ARAIM_USRMASK_GLO;
    abv_mask_gal = measures_all(idusrgal,9) >= ARAIM_USRMASK_GAL;
    abv_mask_bdu = measures_all(idusrbdu,9) >= ARAIM_USRMASK_BDU;
    
    abv_mask = [idusrgps(abv_mask_gps);idusrglo(abv_mask_glo);idusrgal(abv_mask_gal);idusrbdu(abv_mask_bdu)];
    
% Dual frequency measurement noise for integrity
    
    el = measures_all(abv_mask,9);
    prn_index = prn(abv_mask);
    indgps = find(prn_index >= MOPS_MIN_GPSPRN & prn_index <= MOPS_MAX_GPSPRN);
    indglo = find(prn_index >= MOPS_MIN_GLOPRN & prn_index <= MOPS_MAX_GLOPRN);
    indgal = find(prn_index >= MOPS_MIN_GALPRN & prn_index <= MOPS_MAX_GALPRN);
    indbdu = find(prn_index >= MOPS_MIN_BDUPRN & prn_index <= MOPS_MAX_BDUPRN);
        
    ura = ones(length(abv_mask),1)*Inf;
    sig2_if = ones(length(abv_mask),1)*Inf;
    bnom_i = ones(length(abv_mask),1)*Inf;
    psat_i = ones(length(abv_mask),1);  
    
    psat_i(indgps) = ARAIM_PSAT_GPS;
    psat_i(indglo) = ARAIM_PSAT_GLO;
    psat_i(indgal) = ARAIM_PSAT_GAL;
    psat_i(indbdu) = ARAIM_PSAT_BDU;
           
    ura(indgps) = ARAIM_URA_GPS;
    ura(indglo) = ARAIM_URA_GLO;
    ura(indgal) = ARAIM_URA_GAL;
    ura(indbdu) = ARAIM_URA_BDU;
           
    bnom_i(indgps) = ARAIM_BIAS_GPS;
    bnom_i(indglo) = ARAIM_BIAS_GLO;
    bnom_i(indgal) = ARAIM_BIAS_GAL;
    bnom_i(indbdu) = ARAIM_BIAS_BDU;
          
    el = deg2rad(el);
    sig_trv = 0.12*1.001./sqrt(0.002001+sin(el).^2);
            
    dual_freq = 0;  
    
    if flag_sig2_if == 1
        if dual_freq == 1
            sig2_if(indgps) = GAMMA_L1_L5*af_cnmp_mops(NaN,el(indgps));
            sig2_if(indgal) = af_cnmp_galileo(el(indgal));
            sig2_if(indglo) = GAMMA_L1_L5*af_cnmp_mops(NaN,el(indglo)); %use GPS cnmp characterization 
            sig2_if(indbdu) = GAMMA_L1_L5*af_cnmp_mops(NaN,el(indbdu)); %use GPS cnmp characterization 
        else    
            r_ecef = est_r_eb_e + est_C_b_e * L_ba_b; 
            [r_geo(1),r_geo(2),r_geo(3)] = ecef2geo(r_ecef(1),r_ecef(2),r_ecef(3));
            % ll_usr = [rad2deg(r_geo(1)),rad2deg(r_geo(2))];
            ll_usr = [r_geo(1),r_geo(2)];
            % calcualte the IPP latitudes and longitudes
            ll_ipp = calc_ll_ipp(ll_usr.*ones(no_GNSS_meas,1), el_s, az_s);
            % convert back to degrees
            ll_ipp = ll_ipp*180/pi;
            mag_lat = ll_ipp(:,1) + 0.064*180*cos((ll_ipp(:,2)/180-1.617)*pi);
                   
            %mid-latitude klobuchar confidence
            sig2_uire(abv_mask) = 20.25*obliquity2(el);  
            %low-latitude klobuchar confidence
            idx = find(abs(mag_lat) < 20);
            if(~isempty(idx))
                sig2_uire(abv_mask(idx)) = 81*obliquity2(el(idx));
            end
            %high-latitude klobuchar confidence
            idx = find(abs(mag_lat) > 55);
            if(~isempty(idx))
                sig2_uire(abv_mask(idx)) = 36*obliquity2(el(idx));
            end            
            sig2_uire = sig2_uire(abv_mask)';
              
            if dual_freq ==0
                sig2_if(indgps) = sig2_uire(indgps) + af_cnmp_mops(NaN,el(indgps));
                sig2_if(indgal) = sig2_uire(indgal) + af_cnmp_galileo(el(indgal))/GAMMA_L1_L5;
                sig2_if(indglo) = sig2_uire(indglo) + af_cnmp_mops(NaN,el(indglo)); %use GPS cnmp characterization 
                sig2_if(indbdu) = sig2_uire(indbdu) + af_cnmp_mops(NaN,el(indbdu));
        
            else
                sig2_if(indgps) = CONST_GAMMA_L1L5^2*sig2_uire(indgps) + af_cnmp_mops(NaN,el(indgps));
                sig2_if(indgal) = CONST_GAMMA_L1L5^2*sig2_uire(indgal) + af_cnmp_galileo(el(indgal))/GAMMA_L1_L5;
                sig2_if(indglo) = CONST_GAMMA_L1L5^2*sig2_uire(indglo) + af_cnmp_mops(NaN,el(indglo)); %use GPS cnmp characterization 
                sig2_if(indbdu) = CONST_GAMMA_L1L5^2*sig2_uire(indbdu) + af_cnmp_mops(NaN,el(indbdu)); %use GPS cnmp characterization 
            end % end dual_freq==0
        
        end % end dual_freq==1      
    else
        sig2_if = ones(length(abv_mask),1)*0;
    end  % end flag_sig2_if == 1
    
    sig2    = ura.^2+ sig_trv.^2+ sig2_if;


% Dual frequency measurement noise for accuracy
           
    ure = ones(length(abv_mask),1)*Inf;
    bcont_i = ones(length(abv_mask),1)*Inf;
           
    ure(indgps) = ARAIM_URE_GPS;
    ure(indglo) = ARAIM_URE_GLO;
    ure(indgal) = ARAIM_URE_GAL;
    ure(indbdu) = ARAIM_URE_BDU;
           
    bcont_i(indgps) = ARAIM_BIAS_CONT_GPS;
    bcont_i(indglo) = ARAIM_BIAS_CONT_GLO;
    bcont_i(indgal) = ARAIM_BIAS_CONT_GAL;
    bcont_i(indbdu) = ARAIM_BIAS_CONT_BDU;
                
    sig2acc = ure.^2+ sig_trv.^2+ sig2_if;
           

    usridx = ones(length(abv_mask),1);
    [vhpl,lamda,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,Z_all,H_all] = my_usr_vhpl_araim_8(usridx, sig2, ...
        prn_index,...
        sig2acc, bnom_i, bcont_i, psat_i,...  
        measures_all,no_GNSS_meas,tor_s,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
        meas_f_ib_b,TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch,lamda,Z_all,H_all,GNSS_epoch);
    
    bad_usr=find(vhpl(:,1) <= 0 | vhpl(:,4) <= 0);
    if(~isempty(bad_usr))
        vhpl(bad_usr,:)=NaN;
    end

  
end








