clc;clear;
%% Intro
% Code created For ansering the first prompt of AE_313 Final Project. 
% By Shreyas Madhvaraju and Steven Hodac

%% PART a

% Constants and Derivations

mu = 398600;

I_hat = [1 0 0];
J_hat = [0 1 0];
K_hat = [0 0 1];

r_vec = [-15865 13312 41357];
v_vec = [-0.9601 -1.1443 0];

r = norm(r_vec);
v = norm(v_vec);
r_hat = r_vec./r;

h_vec = cross(r_vec,v_vec);
h = norm(h_vec);
h_hat = h_vec./h;

n_hat = cross(K_hat,h_vec)./norm(cross(K_hat,h_vec));

eps = 0.5*v^2-mu/r;

% Conversion to OE at t_1

a = -mu/(2*eps); %semi major axis of orbit
fprintf('Semi major axis of orbit = %.3f (km)\n', a);

e_vec = (1/mu) * cross(v_vec,h_vec) - r_hat;
e = norm(e_vec); %Eccentricty of orbit
fprintf('Eccentricity of orbit = %.3f\n', e);
e_hat = e_vec/e;

i = acos(dot(h_hat,K_hat)); %Inclination of orbit in radians converted to degrees
fprintf('Inclination of orbit = %.3f (deg)\n', i*180/pi);

Omega = atan2(n_hat(2),n_hat(1)); %Longitude of the ascending node
fprintf('Longitude of the ascending node = %.3f (deg)\n', Omega*180/pi);

if dot(h_hat,K_hat) >= 0 %Argument of Periapsis
    w = acos(dot(n_hat,e_hat)); 
    fprintf('Argument of Periapsis = %.3f (deg)\n', w*180/pi);
else 
    w = 2*pi - acos(dot(n_hat,e_hat));
    fprintf('Argument of Periapsis = %.3f (deg)\n', w*180/pi);
end

if dot(r_vec,v_vec) >= 0 %True Anomaly
    f1 = acos(dot(e_hat,r_hat));
    fprintf('Initial True Anomaly = %.3f (deg)\n', f1*180/pi);
else 
    f1 = 2*pi - acos(dot(e_hat,r_hat));
    fprintf('Initial True Anomaly = %.3f (deg)\n', f1*180/pi);
end
pause 
%% PART b

% Newton's method for elipse
if e > 0 && e < 1 
    n = sqrt(mu/a^3);
    E = 2*atan(tan(f1/2) * sqrt((1-e)/(1+e)));
    M_e = E - e*sin(E);
    t_1 = M_e/n;
    fprintf('Time passed since periapsis is %.3f (s)\n' , t_1);
end

% Newton's method for hyperbola
if e > 1
    n = sqrt(-mu/a^3);
    F = 2*atanh(tan(f1/2)*sqrt((e-1)/(e+1)));
    M_h = e*sinh(F)-F;
    t_1 = M_h/n;
    fprintf('Time passed since periapsis is %.3f (s)\n' , t_1);
end
pause
%% PART c

fprintf('Specific energy of orbit is %.3f (km^2/s^2)\n' , eps);
fprintf('Angular momentum vector is [%s] (km^2/s)\n', join(string(h_vec), ','));
fprintf('Angular momentum magnitude is %.3f (km^2/s)\n' , h);
pause

