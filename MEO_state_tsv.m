function [positionX,positionY,positionZ,Vx,Vy,Vz,sat_Clk_error] = MEO_state_tsv(tsv,EPH)

%%%%%%%%%%%%%%%%
radv = 7.2921151467e-5;
GM = 3.986005e14;
F = -4.442807633e-10;
toc = EPH(3,:);
Tsv = tsv;
af0 = EPH(4,:);
af1 = EPH(5,:);
af2 = EPH(6,:);
TGD1 = EPH(29,:);
sqrtA = EPH(14,:);
toe = EPH(15,:);
deltan = EPH(9,:);
M0 = EPH(10,:);
e = EPH(12,:);
i0 = EPH(19,:);
IDOT = EPH(23,:);
omega = EPH(21,:);
Cuc = EPH(11,:);
Cus = EPH(13,:);
Crc = EPH(20,:);
Crs = EPH(8,:);
Cic = EPH(16,:);
Cis = EPH(18,:);
omega0 = EPH(17,:);
omegaDot = EPH(22,:);

delt = Tsv - toc;
sat_Clk_error = af2*delt^2 + af1*delt + af0 - TGD1*1e-9;
Tsv = Tsv - sat_Clk_error;

tk = Tsv - toe;

a = sqrtA^2;
n = sqrt(GM/(a^3)) + deltan;
M = M0 + n*tk;
M = mod(M,2*pi);

E0 = M;
while 1
    Ei = M + e*sin(E0);
    delt_E = Ei - E0;
    delt_E = rem(delt_E,2*pi);
    if abs(delt_E) < 1e-12
        E = Ei;
        break;
    else
        E0 = Ei;
    end
end

nu = atan2(sqrt(1-e^2)*sin(E),cos(E)-e);
phi = nu + omega;
phi = rem(phi,2*pi);

delt_u = Cuc*cos(2*phi) + Cus*sin(2*phi);
delt_r = Crc*cos(2*phi) + Crs*sin(2*phi);
delt_i = Cic*cos(2*phi) + Cis*sin(2*phi);
u = phi + delt_u;
r = a*(1-e*cos(E))+delt_r;
i = i0+IDOT*tk+delt_i;

OMEGA = omega0+(omegaDot-radv)*tk-radv*toe;
OMEGA = mod(OMEGA,2*pi);

positionX = cos(u)*r*cos(OMEGA) - sin(u)*r*cos(i)*sin(OMEGA);
positionY = cos(u)*r*sin(OMEGA) + sin(u)*r*cos(i)*cos(OMEGA);
positionZ = sin(u)*r*sin(i);
 

E_dot = n/(1-e*cos(E));
Vk_dot = (sqrt(1-e^2)*E_dot)/(1-e*cos(E));
 
Uk_dot = 2*Vk_dot*(Cus*cos(2*phi)-Cuc*sin(2*phi));
Rk_dot = 2*Vk_dot*(Crs*cos(2*phi)-Crc*sin(2*phi));
Ik_dot = 2*Vk_dot*(Cis*cos(2*phi)-Cic*sin(2*phi));
 
U_dot = Vk_dot + Uk_dot;
R_dot = a*e*E_dot*sin(E) + Rk_dot;
I_dot = IDOT + Ik_dot;
 
OMEGA_dot = omegaDot - radv;
 
x = R_dot*cos(u) - r*U_dot*sin(u);
y = R_dot*sin(u) + r*U_dot*cos(u);
 
Vx = -positionY*OMEGA_dot - (y*cos(i)-positionZ*I_dot)*sin(OMEGA) + x*cos(OMEGA);
Vy = positionX*OMEGA_dot + (y*cos(i)-positionZ*I_dot)*cos(OMEGA) + x*sin(OMEGA);
Vz = y*sin(i) + y*I_dot*cos(i);
 
dtr = F*e*sqrtA*sin(E);
sat_Clk_error = af2*delt^2 + af1*delt + af0 - TGD1*1e-9 + dtr;
 
 
end
 
