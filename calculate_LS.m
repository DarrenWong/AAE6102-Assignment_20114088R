clc;
clear all;
close all;
format long g
% define the constant
% speed of light
c = 299792458.0;
% earth’s rotation rate (r/s)
Wedot= 7.2921151467e-5;
% WGS 84 value of earth's universal gravitation constant: (m^3/s^2);
mu= 3.986005e+14;
% Relativistic correction term constant
F_kep= (-4.442807633e-10);
      
XR0 = [-2694685.473 -4293642.366 3857878.924];
reference_result = [-2700400 -4292560 3855270];
% load the ephemeris and observation data
load('Data/eph.dat')
load('Data/rcvr.dat')

user_pose_init = XR0;
receiver_time = rcvr(1,1);
pr_R   = zeros(32,1);
elR   = zeros(32,1);

%read obs data
for i = 1 :32
    k = find(rcvr(:,2)==i);
    if(isempty(k)) 
       continue;
    end
    pr_R(i) = rcvr(k,3) ;
end
index = find(pr_R~=0);

% init para for the iterations
n_iter=0;
dXR = Inf;
XR_data = zeros(10,3); 
dtR = 0;

% sat position
[XS,dtS] = sat_position(receiver_time,pr_R, eph);

%condition for LS estimation, dXR > 10^-4 or iteration larger 10
while(dXR >10^-4 || n_iter>10)
    
    n_iter = n_iter + 1;

    % save initial 
    XR_data(n_iter,:) = XR0;
    
    %tropo_error_correction
    [wlat, wlon, walt] = Wgsxyz2lla(XR0);
    for jdx = 1 : length(XS)
        if(~all(XS(jdx)))
            continue;
        end
        target_pos_xyz = [XS(jdx,:)];    
        target_pos_enu = xyz2enu(target_pos_xyz,XR0);
        elR(jdx) = ( 180/pi*atan(target_pos_enu(3)/norm(target_pos_enu(1:2))));
    end
    err_tropo_R = tropo_error_correction(elR(elR ~= 0), walt);
   
    % compute least square
    [XR, dtR, dXR] = LS_code(XR0, XS(index,:), pr_R(index),dtS(index),dtR, err_tropo_R);
    XR0=XR;    
    fprintf("iteration %d, pose error %d user clock bias %d, dXR %d\n",n_iter, norm(XR-reference_result),dtR, dXR);
end  

fprintf("pose estimation result\n");
fprintf("total iteration %d, pose error %d user clock bias%d\n",n_iter, norm(XR-reference_result),dtR);

% receive position ECEF to LLA 
[wlat, wlon, walt] = Wgsxyz2lla(XR);
% inti position ECEF to LLA 

% init_pos_enu =Wgsxyz2enu(user_pose_init',wlat, wlon, walt);
% reference_pos_enu =Wgsxyz2enu(reference_result',wlat, wlon, walt);
% 
% figure(1)
% geoplot(wlat,wlon,'.','Color','b','MarkerSize',10)
% geolimits([wlat-0.1 wlat+0.1],[wlon-0.1 wlon+0.1])
% geobasemap streets
% hold on
% 
% figure(2)
% % x = [init_pos_enu(1),0,reference_pos_enu(1)]
% % y = [init_pos_enu(2),0,reference_pos_enu(2)]
% 
% plot(init_pos_enu(1),init_pos_enu(2),'.','MarkerSize',30)
% ylim([-20 4000])
% xlim([-20 4000])
% hold on
% plot(0,0,'.','MarkerSize',30)
% hold on
% plot(reference_pos_enu(1),reference_pos_enu(2),'^','MarkerSize',30)


% plot(init_pos_enu(1),init_pos_enu(2),'*r')
% ylim([-20 4000])
% xlim([-20 4000])
% hold on
% plot(0,0,'*b')



function [XS, dtS] = sat_position(t, pseudorange, eph)
% INPUT:
%   t =  reception time
%   pseudorange = observed pseudoranges
%   Eph = ephemerides matrix
% OUTPUT:
%   XS      = satellite position at transmission time in ECEF(time_rx) (X,Y,Z)
%   dtS     = satellite clock error (vector)
%
% DESCRIPTION:
%   Computation of the satellite positions based on reception time.

  XS_tx   = zeros(32,3);
  XS   = zeros(32,3);
  dtS   = zeros(32,1);
  
  c = 299792458.0;
  % earth’s rotation rate (r/s)
  Wedot= 7.2921151467e-5;
  % WGS 84 value of earth's universal gravitation constant: (m^3/s^2);
  mu= 3.986005e+14;
  % Relativistic correction term constant
  F_kep= (-4.442807633e-10);
  nsat = 32;
  for i = 1 :nsat
    k = find(eph(:,2)==i);
    if(isempty(k)) 
       continue;
    end
    Eph = eph(k,:);
    rcvr_tow  = Eph(1);
    svid      = Eph(2);
    toc       = Eph(3);
    toe       = Eph(4);
    af0       = Eph(5);
    af1       = Eph(6);
    af2       = Eph(7);
    ura       = Eph(8);
    e         = Eph(9);
    sqrta     = Eph(10);
    dn        = Eph(11);
    m0        = Eph(12);
    omg       = Eph(13);
    omg0      = Eph(14);
    i0        = Eph(15);
    odot      = Eph(16);
    idot      = Eph(17);
    cus       = Eph(18);
    cus       = Eph(18);
    cuc       = Eph(19);
    cis       = Eph(20);
    cic       = Eph(21);
    crs       = Eph(22);
    crc       = Eph(23);
    iod       = Eph(24);
    % transmission time by satellite clock and pseudorange
    %rcvr(find(rcvr(:,2)==i),3);
    time_tx_raw =t - pseudorange(i)/c;
    Ek = ecc_anomaly(time_tx_raw, Eph);
    relativistic_error_correction = F_kep * e* sqrta *sin(Ek);
    
    deltaS = af0 + af1*(time_tx_raw-toc) + af2*(time_tx_raw-toc)^2;
    deltaS = deltaS + relativistic_error_correction;
    % corrected sat clock
    time_tx_corr = time_tx_raw - deltaS;
    % satellite clock error
    deltaS = af0 + af1*(time_tx_corr-toc) + af2*(time_tx_corr-toc)^2;
    deltaS = deltaS + relativistic_error_correction;
    dtS(i,:) =  deltaS;
    time_tx = time_tx_raw - dtS(i,:);
   %semi-major axis
    A = sqrta*sqrta; 
    % Calculate the satellite position ICD
    tk = time_tx - toe;
    %eccentric anomaly
    Ek = ecc_anomaly(time_tx, Eph);
    % true anomaly
    vk = atan2(sqrt(1-e^2)*sin(Ek),cos(Ek)-e);
    % argument of Latitude
    phik = vk + omg;
    phik = rem(phik,2*pi); 
    % second harmonic perturbations 
    %  corrected argument of latitude
    uk = phik + cuc*cos(2*phik) + cus*sin(2*phik);
    %  corrected radius 
    rk = A*(1 - e*cos(Ek)) + crc*cos(2*phik) + crs*sin(2*phik);  
    % corrected inclination
    ik = i0 + cic*cos(2*phik) + cis*sin(2*phik)+ idot*tk ;  
    % position in orbital plane
    xk_prime = rk * cos(uk);
    yk_prime = rk * sin(uk);
    %corrected longitude of ascending node
    omegak = omg0 + (odot - Wedot)*tk - Wedot*toe;
    omegak = rem(omegak + 2*pi, 2*pi);
      
    %satellite Earth-fixed coordinates (X,Y,Z)
    xk = xk_prime*cos(omegak) - yk_prime*cos(ik)*sin(omegak);
    yk = xk_prime*sin(omegak) + yk_prime*cos(ik)*cos(omegak);
    zk = yk_prime*sin(ik);

    XS_tx(i,:) = [xk yk zk];

    %column data
    Xsat = transpose(XS_tx(i,:));
   % time_tx=0;
    %computation of ECEF satellite position at time_rx
    traveltime = t -time_tx;  

    %find the rotation angle
    omegatau = Wedot* traveltime;
    %build a rotation matrix
    Rot = [ cos(omegatau)    sin(omegatau)   0;
          -sin(omegatau)    cos(omegatau)   0;
           0                0               1];
    %Correct satellite position according to Earth rotation during signal travel time.
    XS(i,:) =Rot*Xsat;
  end
    
end

function [Ek, n] = ecc_anomaly(time, Eph)
% INPUT:
%   time = GPS time
%   Eph = ephemerides matrix
% OUTPUT:
%   Ek = eccentric anomaly
%   n = corrected mean motion [rad/sec]
%
% DESCRIPTION:
%   Computation of the eccentric anomaly.

    mu= 3.986005e+14;
    %get ephemerides
    m0        = Eph(12);
    sqrta     = Eph(10);
    dn        = Eph(11);
    e         = Eph(9);
    toe       = Eph(4);
    
    tk = time - toe;
    % semi-major axis
    A= sqrta * sqrta;
    % computed mean motion(rad/sec)
    n0 = sqrt(mu/A^3);
    % corrected mean motion
    n = n0 + dn;
    % mean anomaly
    Mk = m0 + n*tk;
    % Eccentric Anomaly by iteration Mk=Ek-esinEk, refer to gogps and satpos.m Kai Borre
    Mk = rem(Mk+2*pi,2*pi); % between 0 and 360 
    Ek = Mk;   %Initial guess of eccentric anomaly
    for idx = 1 : 10
       Ek_old = Ek;
       Ek = Mk+e*sin(Ek);
       dEk = rem(Ek-Ek_old,2*pi);
       if abs(dEk) < 1.e-12
          break
       end
    end

    if (idx == 11)
        fprintf('WARNING: Eccentric anomaly does not converge.\n');
    end
    Ek = rem(Ek+2*pi,2*pi);  % between 0 and 360 
end 

function [XR, dtR, sigma_XR] = LS_code(XR_approx, XS, pr_R,dtS, dtR, err_tropo_RS)
% INPUT:
%   el = satellite elevation
%   h  = receiver ellipsoidal height
%
% OUTPUT:
%   XR = tropospheric error correction
% DESCRIPTION:
%   least squares to estimate receiver location 

    v_light = 299792458.0;

    n = length(pr_R);

    XR_mat = XR_approx(ones(n,1),:);

    distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));

    %known term vector
    b = distR_approx + v_light*(dtR - dtS) + err_tropo_RS; 

    %observation vector
    y0 = pr_R;
  
   %design matrix H
    H = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
         (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
         (XR_approx(3) - XS(:,3)) ./ distR_approx, ... %column for Z coordinate
          ones(n,1)];        %column for receiver clock delay (multiplied by c) 

    dx =  (H'*H)^(-1)  * H' * (y0-b);
    XR =  XR_approx + dx(1:3)';
    dtR = dtR + dx(4)/v_light;
    sigma_XR = norm(dx(1:3));
end

function [corr] = tropo_error_correction(el, h)
% SYNTAX:
%   [corr] = tropo_error_correction(el, h);
%
% INPUT:
%   el = satellite elevation
%   h  = receiver ellipsoidal height
%
% OUTPUT:
%   corr = tropospheric error correction
%
% DESCRIPTION:
%   Computation of the pseudorange correction due to tropospheric refraction.
%   Saastamoinen algorithm.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%Saastamoinen model requires positive ellipsoidal height
    h(h < 0) = 0;
	%conversion to radians
	el = abs(el) * pi/180;

	    %Standard atmosphere - Berg, 1948 (Bernese)
	%pressure [mbar]
	Pr = 1013.25;
	%temperature [K]
	Tr = 291.15;
	%numerical constants for the algorithm [-] [m] [mbar]
	Hr = 50.0;

	P = Pr * (1-0.0000226*h).^5.225;
	T = Tr - 0.0065*h;
	H = Hr * exp(-0.0006396*h);

	%----------------------------------------------------------------------

	%linear interpolation
	h_a = [0; 500; 1000; 1500; 2000; 2500; 3000; 4000; 5000];
	B_a = [1.156; 1.079; 1.006; 0.938; 0.874; 0.813; 0.757; 0.654; 0.563];

	t = zeros(length(T),1);
	B = zeros(length(T),1);

	for i = 1 : length(T)

	    d = h_a - h(i);
	    [dmin, j] = min(abs(d));
	    if (d(j) > 0)
		index = [j-1; j];
	    else
		index = [j; j+1];
	    end

	    t(i) = (h(i) - h_a(index(1))) ./ (h_a(index(2)) - h_a(index(1)));
	    B(i) = (1-t(i))*B_a(index(1)) + t(i)*B_a(index(2));
    end
	%----------------------------------------------------------------------
	e = 0.01 * H .* exp(-37.2465 + 0.213166*T - 0.000256908*T.^2);

	%tropospheric error
	corr = ((0.002277 ./ sin(el)) .* (P - (B ./ (tan(el)).^2)) + (0.002277 ./ sin(el)) .* (1255./T + 0.05) .* e);
end

