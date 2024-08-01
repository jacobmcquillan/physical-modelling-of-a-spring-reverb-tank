%% SET MATERIAL & GEOMETRY OF SPRING TANK %%
% JACOB MCQUILLAN
% SARC, QUEEN'S UNIVERSITY BELFAST

function[] = reverbTankParams()
%% SEGMENT 5: HELICAL SPRING %%
% key
% L5: unwound wire length [metres] = circumference of coil x no. of coils 
% r5: wire radius [metres]
% Rc5: coil radius [metres]
% alp5: helix angle [rad] = arctan(pitch/coil diameter),
% pitch = free length / no. of coils

% olson (Parker, 2009) - CHANGE MAGNETIC BEAD DIMENSIONS TOO
Lf = 0.065; coils = 148;
r5 = 0.175e-3; Rc5 = 2.7e-3;
E5 = 2e11; % Young's modulus [Pa]
nu5 = 0.3; % Poisson's ratio
rho5 = 7800; % material density [kg/m^3]

pitch = Lf/coils;
L5 = 2*pi*Rc5*coils;
alp5 = atan(pitch/(2*pi*Rc5));

% % spring 2 in paper
% E5 = 2.27e11; % Young's modulus [Pa]
% nu5 = 0.3; % Poisson's ratio
% rho5 = 7800; % material density [kg/m^3]
% r5 = 0.15e-3;
% Rc5 = 2.3e-3 - r5;
% coils = 352;
% Lf = 0.163;
% 
% pitch = Lf/coils;
% L5 = 2*pi*Rc5*coils;
% alp5 = atan(pitch/(2*pi*Rc5));

% % spring 1 in paper
% E5 = 2.32e11; % Young's modulus [Pa]
% nu5 = 0.3; % Poisson's ratio
% rho5 = 7800; % material density [kg/m^3]
% 
% r5 = 0.15e-3;
% Rc5 = 2.3e-3 - r5;
% coils = 297;
% Lf = 0.163;
% 
% pitch = Lf/coils;
% L5 = 2*pi*Rc5*coils;
% alp5 = atan(pitch/(2*pi*Rc5));

% derived parameters
G5 = E5 / (2*(1+nu5)); % Shear modulus [Pa]
I5 = (pi*(r5^4)) / 4; % 2nd moment of area
CSA5 = pi*(r5^2); % cross-sectional area

kap5 = ((cos(alp5))^2) / Rc5; % helix curvature
mu5 = tan(alp5); % kap*mu = tau is helix tortuosity


%% SEGMENT 1: WIRE AT LEFT ENDPOINT %%
% berrylium copper (Young, 1961)
E1 = 1.7e11; nu1 = 0.3; rho1 = 8000;
r1 = 0.1e-3;
L1 = 20e-3;

G1 = E1 / (2*(1+nu1)); I1 = (pi*(r1^4)) / 4; CSA1 = pi*(r1^2);

gam1u = (E1*I1) / (rho1*CSA1);
gam1v = gam1u;
gam1w = E1 / rho1;
gam1thet = G1 / rho1;


%% SEGMENT 9: WIRE AT RIGHT ENDPOINT %%
% same as segment 1
E9 = E1; nu9 = nu1; rho9 = rho1;
r9 = r1; L9 = L1;

G9 = E9 / (2*(1+nu9)); I9 = (pi*(r9^4)) / 4; CSA9 = pi*(r9^2);

gam9u = (E9*I9) / (rho9*CSA9);
gam9v = gam9u;
gam9w = E9 / rho9;
gam9thet = G9 / rho9;


%% SEGMENT 2: LEFT MAGNETIC BEAD %%
% ceramic magnetic material (Meinema, 1961)
E2 = 1.5e8; rho2 = 5000; nu2 = 0.3;

% spring tank in paper
% r2 = 0.9e-3;
% L2 = 4e-3;

% olson (Parker, 2009)
L2 = 5e-3;
r2 = 1.05e-3;

G2 = E2 / (2*(1+nu2)); I2 = (pi*(r2^4)) / 4; CSA2 = pi*(r2^2);

gam2u = (E2*I2) / (rho2*CSA2);
gam2v = gam2u;
gam2w = E2 / rho2;
gam2thet = G2 / rho2;


%% SEGMENT 8: RIGHT MAGNETIC BEAD
% same as segment 2
E8 = E2; nu8 = nu2; rho8 = rho2;
r8 = r2; L8 = L2;

G8 = E8 / (2*(1+nu8)); I8 = (pi*(r8^4)) / 4; CSA8 = pi*(r8^2);

gam8u = (E8*I8) / (rho8*CSA8);
gam8v = gam8u;
gam8w = E8 / rho8;
gam8thet = G8 / rho8;


%% SEGMENT 3: WIRE AFTER LEFT BEAD %%
E3 = E5; nu3 = nu5; rho3 = rho5;
r3 = r5; L3 = 5e-3;

G3 = E3 / (2*(1+nu3)); I3 = (pi*(r3^4)) / 4; CSA3 = pi*(r3^2);

gam3u = (E3*I3) / (rho3*CSA3);
gam3v = gam3u;
gam3w = E3 / rho3;
gam3thet = G3 / rho3;


%% SEGMENT 7: WIRE BEFORE RIGHT BEAD %%
% same as segment 3
E7 = E3; nu7 = nu3; rho7 = rho3;
r7 = r3; L7 = L3;

G7 = E7 / (2*(1+nu7)); I7 = (pi*(r7^4)) / 4; CSA7 = pi*(r7^2);

gam7u = (E7*I7) / (rho7*CSA7);
gam7v = gam7u;
gam7w = E7 / rho7;
gam7thet = G7 / rho7;


%% SEGMENT 4: WIRE BEFORE SPRING %%
E4 = E5; nu4 = nu5; rho4 = rho5;
r4 = r5; L4 = Rc5;

G4 = E4 / (2*(1+nu4)); I4 = (pi*(r4^4)) / 4; CSA4 = pi*(r4^2);

gam4u = (E4*I4) / (rho4*CSA4);
gam4v = gam4u;
gam4w = E4 / rho4;
gam4thet = G4 / rho4;


%% SEGMENT 6: WIRE AFTER SPRING %%
% same as segment 4
E6 = E5; nu6 = nu5; rho6 = rho5;
r6 = r5; L6 = L4;

G6 = E6 / (2*(1+nu6)); I6 = (pi*(r6^4)) / 4; CSA6 = pi*(r6^2);

gam6u = (E6*I6) / (rho6*CSA6);
gam6v = gam6u;
gam6w = E6 / rho6;
gam6thet = G6 / rho6;


%% CONSOLIDATE EVERYTHING TO KEEP TIDY WORKSPACE %%
E = [E1 E2 E3 E4 E5 E6 E7 E8 E9];
nu = [nu1 nu2 nu3 nu4 nu5 nu6 nu7 nu8 nu9];
rho = [rho1 rho2 rho3 rho4 rho5 rho6 rho7 rho8 rho9];

L = [L1 L2 L3 L4 L5 L6 L7 L8 L9];
r = [r1 r2 r3 r4 r5 r6 r7 r8 r9];

G = [G1 G2 G3 G4 G5 G6 G7 G8 G9];
I = [I1 I2 I3 I4 I5 I6 I7 I8 I9];
CSA = [CSA1 CSA2 CSA3 CSA4 CSA5 CSA6 CSA7 CSA8 CSA9];

gamu = [gam1u gam2u gam3u gam4u 0 gam6u gam7u gam8u gam9u];
gamv = [gam1v gam2v gam3v gam4v 0 gam6v gam7v gam8v gam9v];
gamw = [gam1w gam2w gam3w gam4w 0 gam6w gam7w gam8w gam9w];
gamthet = [gam1thet gam2thet gam3thet gam4thet 0 gam6thet gam7thet gam8thet gam9thet];


%% UPDATE VALUES %%
save('SpringReverbParameters.mat','E','nu','rho','L','r','G','I','CSA','gamu','gamv','gamw','gamthet','Rc5','alp5','kap5','mu5','coils','Lf');