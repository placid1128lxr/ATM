function varargout=orbital_propagator_BDS0(eph,tx_raw)

%number of satellites
k=size(eph,2);
positionX = zeros(k,1);
positionY = zeros(k,1);
positionZ = zeros(k,1);
sat_Clk_error = zeros(k,1);

idx_GEO = [1:5,59:62];
prn = eph(1,:);
for i = 1:k
    if(ismember(prn(i),idx_GEO))
        [positionX(i),positionY(i),positionZ(i),sat_Clk_error(i)] = GEO_state_tsv(tx_raw(i),eph(:,i));
    else
        [positionX(i),positionY(i),positionZ(i),~,~,~,sat_Clk_error(i)] = MEO_state_tsv(tx_raw(i),eph(:,i));
    end
end

position_ecef = [positionX,positionY,positionZ]';
dt_sv = sat_Clk_error;

%outputs
varargout{1}=position_ecef;
varargout{2}=dt_sv;

    



