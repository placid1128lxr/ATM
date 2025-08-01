function [positionX,positionY,positionZ,sat_Clk_error] = GEO_state_tsv(tsv,EPH)

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

OMEGA = omega0+omegaDot*tk-radv*toe;
OMEGA = mod(OMEGA,2*pi);

X = cos(u)*r*cos(OMEGA) - sin(u)*r*cos(i)*sin(OMEGA);
Y = cos(u)*r*sin(OMEGA) + sin(u)*r*cos(i)*cos(OMEGA);
Z = sin(u)*r*sin(i);

f = -5/180*pi;
p = radv*tk;
Rx = [1 0 0; 0 cos(f) sin(f); 0 -sin(f) cos(f)];
Rz = [cos(p) sin(p) 0; -sin(p) cos(p) 0; 0 0 1];
position = Rz*Rx*[X;Y;Z];
positionX = position(1);
positionY = position(2);
positionZ = position(3);

dtr = F*e*sqrtA*sin(E);
sat_Clk_error = af2*delt^2 + af1*delt + af0 - TGD1*1e-9 + dtr;
 
end
 
