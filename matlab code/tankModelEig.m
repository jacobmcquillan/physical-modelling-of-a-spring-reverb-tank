%% FORM SPRING REVERB TANK MODEL EIGENVALUE PROBLEM %%
% JACOB MCQUILLAN
% SARC, QUEEN'S UNIVERSITY BELFAST

function[P,Q,s,Zw,Zr] = tankModelEig(f_BCs,N,eta,sigA,sigB,sigC)
%% LOAD MATERIAL PARAMETERS %%
reverbTankParams(); % change model parameters in reverbTankParams.m
load SpringReverbParameters


%% KELVIN-VOIGT DAMPING %%
for n=1:9 % depends on Young's modulus of each segment
    psi(n) = eta/E(n);
end


%% ANGLED CONNECTION %%
A = sin(alp5); B = cos(alp5);


%% SETUP FOR SPATIAL DISCRETISATION %%
% parse grid points
N1 = N(1); N9 = N(9); % wires between endpoints and beads
N2 = N(2); N8 = N(8); % magnetic beads
N3 = N(3); N7 = N(7); % wires after beads
N4 = N(4); N6 = N(6); % wires at springs
N5 = N(5); % helical spring

% convenient indexing
M1 = N1;
M2 = M1 + N2;
M3 = M2 + N3;
M4 = M3 + N4;
M5 = M4 + N5;
M6 = M5 + N6;
M7 = M6 + N7;
M8 = M7 + N8;
M9 = M8 + N9;

% spatial domains
dom1 = [0 L(1)]; dom9 = [0 L(9)];
dom2 = [0 L(2)]; dom8 = [0 L(8)];
dom3 = [0 L(3)]; dom7 = [0 L(7)];
dom4 = [0 L(4)]; dom6 = [0 L(6)];
dom5 = [0 L(5)];

% Gauss-Chebyshev-Lobatto grid points for each segment
s1 = chebpts(N1,dom1); s9 = chebpts(N9,dom9);
s2 = chebpts(N2,dom2); s8 = chebpts(N8,dom8);
s3 = chebpts(N3,dom3); s7 = chebpts(N7,dom7);
s4 = chebpts(N4,dom4); s6 = chebpts(N6,dom6);
s5 = chebpts(N5,dom5);

% concatenate segment grids
s = [s1; s2; s3; s4; s5; s6; s7; s8; s9];


%% BUILD CHEBYSHEV DIFFERENTIATION MATRICES %%
% notation: Dxy - yth order derivative for xth segment

% segment 1
D11 = diffmat(N1,1,dom1); D12 = diffmat(N1,2,dom1);
D13 = diffmat(N1,3,dom1); D14 = diffmat(N1,4,dom1);

% segment 9
D91 = diffmat(N9,1,dom9); D92 = diffmat(N9,2,dom9);
D93 = diffmat(N9,3,dom9); D94 = diffmat(N9,4,dom9);

% segment 2
D21 = diffmat(N2,1,dom2); D22 = diffmat(N2,2,dom2);
D23 = diffmat(N2,3,dom2); D24 = diffmat(N2,4,dom2);

% segment 8
D81 = diffmat(N8,1,dom8); D82 = diffmat(N8,2,dom8);
D83 = diffmat(N8,3,dom8); D84 = diffmat(N8,4,dom8);

% segment 3
D31 = diffmat(N3,1,dom3); D32 = diffmat(N3,2,dom3);
D33 = diffmat(N3,3,dom3); D34 = diffmat(N3,4,dom3);

% segment 7
D71 = diffmat(N7,1,dom7); D72 = diffmat(N7,2,dom7);
D73 = diffmat(N7,3,dom7); D74 = diffmat(N7,4,dom7);

% segment 4
D41 = diffmat(N4,1,dom4); D42 = diffmat(N4,2,dom4);
D43 = diffmat(N4,3,dom4); D44 = diffmat(N4,4,dom4);

% segment 6
D61 = diffmat(N6,1,dom6); D62 = diffmat(N6,2,dom6);
D63 = diffmat(N6,3,dom6); D64 = diffmat(N6,4,dom6);

% segment 5
D51 = diffmat(N5,1,dom5); D52 = diffmat(N5,2,dom5);
D53 = diffmat(N5,3,dom5); D54 = diffmat(N5,4,dom5);


%% BEAM SYSTEM MATRICES %%
% segments 1 & 9
Gam1u = -gamu(1) * D14; Gam1v = -gamv(1) * D14;
Gam1w = gamw(1) * D12; Gam1thet = gamthet(1) * D12;

Gam9u = -gamu(9) * D94; Gam9v = -gamv(9) * D94;
Gam9w = gamw(9) * D92; Gam9thet = gamthet(9) * D92;

% segments 2 & 8
Gam2u = -gamu(2) * D24; Gam2v = -gamv(2) * D24;
Gam2w = gamw(2) * D22; Gam2thet = gamthet(2) * D22;

Gam8u = -gamu(8) * D84; Gam8v = -gamv(8) * D84;
Gam8w = gamw(8) * D82; Gam8thet = gamthet(8) * D82;

% segments 3 & 7
Gam3u = -gamu(3) * D34; Gam3v = -gamv(3) * D34;
Gam3w = gamw(3) * D32; Gam3thet = gamthet(3) * D32;

Gam7u = -gamu(7) * D74; Gam7v = -gamv(7) * D74;
Gam7w = gamw(7) * D72; Gam7thet = gamthet(7) * D72;

% segments 4 & 6
Gam4u = -gamu(4) * D44; Gam4v = -gamv(4) * D44;
Gam4w = gamw(4) * D42; Gam4thet = gamthet(4) * D42;

Gam6u = -gamu(6) * D64; Gam6v = -gamv(6) * D64;
Gam6w = gamw(6) * D62; Gam6thet = gamthet(6) * D62;


%% SPRING SYSTEM MATRIX %%
Imat = eye(N5); % identity matrix

% all of dimensions (N5 x N5)
Gam5uu = ( -E(5)*I(5)*D54 + 6*E(5)*I(5)*D52*kap5^2*mu5^2 + Imat*(-E(5)*I(5)*kap5^4*mu5^4 - 2*G(5)*I(5)*kap5^4*mu5^2 - CSA(5)*E(5)*kap5^2) ) / (CSA(5)*rho(5));
Gam5vu = ( 4*D51*E(5)*I(5)*kap5^3*mu5^3 + 2*G(5)*D51*I(5)*kap5^3*mu5 - 4*D53*E(5)*I(5)*kap5*mu5 ) / (CSA(5)*rho(5));
Gam5wu = ( E(5)*I(5)*D53*kap5 - E(5)*D51*kap5*CSA(5) - 3*E(5)*I(5)*D51*kap5^3*mu5^2 ) / (CSA(5)*rho(5));
Gam5thetu = ( E(5)*D51*kap5^2*mu5 + G(5)*D51*kap5^2*mu5 ) / rho(5);

Gam5uv = ( 4*E(5)*I(5)*D53*kap5*mu5 - 4*E(5)*I(5)*D51*kap5^3*mu5^3 - 2*G(5)*I(5)*D51*kap5^3*mu5 ) / (CSA(5)*rho(5));
Gam5vv = ( -E(5)*I(5)*D54 + 6*E(5)*I(5)*D52*kap5^2*mu5^2 + 2*G(5)*I(5)*D52*kap5^2 - Imat*(E(5)*I(5)*kap5^4*mu5^4) ) / (CSA(5)*rho(5));
Gam5wv = ( Imat*(E(5)*I(5)*kap5^4*mu5^3) - 3*D52*E(5)*I(5)*kap5^2*mu5 ) / (CSA(5)*rho(5));
Gam5thetv = ( E(5)*D52*kap5 - Imat*(E(5)*kap5^3*mu5^2) + 2*G(5)*D52*kap5 ) / (2*rho(5));

Gam5uw = ( CSA(5)*E(5)*D51*kap5 - E(5)*I(5)*D53*kap5 + 3*E(5)*I(5)*D51*kap5^3*mu5^2 ) / (CSA(5)*rho(5));
Gam5vw = ( Imat*(E(5)*I(5)*kap5^4*mu5^3) - 3*D52*E(5)*I(5)*kap5^2*mu5 ) / (CSA(5)*rho(5));
Gam5ww = ( E(5)*I(5)*D52*kap5^2 + CSA(5)*E(5)*D52 - Imat*(mu5^2*E(5)*I(5)*kap5^4) ) / (CSA(5)*rho(5));
Gam5thetw = ( Imat*(E(5)*kap5^3*mu5) ) / (2*rho(5));

Gam5uthet = ( -2*G(5)*I(5)*D51*kap5^2*mu5 - 2*E(5)*I(5)*D51*kap5^2*mu5 ) / (CSA(5)*rho(5));
Gam5vthet = ( E(5)*I(5)*D52*kap5 + 2*G(5)*I(5)*D52*kap5 - Imat*(E(5)*I(5)*kap5^3*mu5^2) ) / (CSA(5)*rho(5));
Gam5wthet = ( Imat*(E(5)*I(5)*kap5^3*mu5) ) / (CSA(5)*rho(5));
Gam5thetthet = ( 2*G(5)*D52 - Imat*(kap5^2*E(5)) ) / (2*rho(5));


%% CONSOLIDATE ALL SYSTEM MATRICES %%
% u polarisation for all segments
GamU = zeros(M9,4*M9);
% beams
GamU(1:M1,1:M1) = Gam1u;
GamU(M1+1:M2,M1+1:M2) = Gam2u;
GamU(M2+1:M3,M2+1:M3) = Gam3u;
GamU(M3+1:M4,M3+1:M4) = Gam4u;
GamU(M5+1:M6,M5+1:M6) = Gam6u;
GamU(M6+1:M7,M6+1:M7) = Gam7u;
GamU(M7+1:M8,M7+1:M8) = Gam8u;
GamU(M8+1:M9,M8+1:M9) = Gam9u;
% spring
GamU(M4+1:M5,M4+1:M5) = Gam5uu;
GamU(M4+1:M5,M9+M4+1:M9+M5) = Gam5uv;
GamU(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = Gam5uw;
GamU(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = Gam5uthet;

% v polarisation for all segments
GamV = zeros(M9,4*M9);
% beams
GamV(1:M1,M9+(1:M1)) = Gam1v;
GamV(M1+1:M2,M9+(M1+1:M2)) = Gam2v;
GamV(M2+1:M3,M9+(M2+1:M3)) = Gam3v;
GamV(M3+1:M4,M9+(M3+1:M4)) = Gam4v;
GamV(M5+1:M6,M9+(M5+1:M6)) = Gam6v;
GamV(M6+1:M7,M9+(M6+1:M7)) = Gam7v;
GamV(M7+1:M8,M9+(M7+1:M8)) = Gam8v;
GamV(M8+1:M9,M9+(M8+1:M9)) = Gam9v;
% spring
GamV(M4+1:M5,M4+1:M5) = Gam5vu;
GamV(M4+1:M5,M9+M4+1:M9+M5) = Gam5vv;
GamV(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = Gam5vw;
GamV(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = Gam5vthet;

% w polarisation for all segments
GamW = zeros(M9,4*M9);
% beams
GamW(1:M1,(2*M9)+(1:M1)) = Gam1w;
GamW(M1+1:M2,(2*M9)+(M1+1:M2)) = Gam2w;
GamW(M2+1:M3,(2*M9)+(M2+1:M3)) = Gam3w;
GamW(M3+1:M4,(2*M9)+(M3+1:M4)) = Gam4w;
GamW(M5+1:M6,(2*M9)+(M5+1:M6)) = Gam6w;
GamW(M6+1:M7,(2*M9)+(M6+1:M7)) = Gam7w;
GamW(M7+1:M8,(2*M9)+(M7+1:M8)) = Gam8w;
GamW(M8+1:M9,(2*M9)+(M8+1:M9)) = Gam9w;
% spring
GamW(M4+1:M5,M4+1:M5) = Gam5wu;
GamW(M4+1:M5,M9+M4+1:M9+M5) = Gam5wv;
GamW(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = Gam5ww;
GamW(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = Gam5wthet;

% thet polarisation for all segments
GamThet = zeros(M9,4*M9);
% beams
GamThet(1:M1,(3*M9)+(1:M1)) = Gam1thet;
GamThet(M1+1:M2,(3*M9)+(M1+1:M2)) = Gam2thet;
GamThet(M2+1:M3,(3*M9)+(M2+1:M3)) = Gam3thet;
GamThet(M3+1:M4,(3*M9)+(M3+1:M4)) = Gam4thet;
GamThet(M5+1:M6,(3*M9)+(M5+1:M6)) = Gam6thet;
GamThet(M6+1:M7,(3*M9)+(M6+1:M7)) = Gam7thet;
GamThet(M7+1:M8,(3*M9)+(M7+1:M8)) = Gam8thet;
GamThet(M8+1:M9,(3*M9)+(M8+1:M9)) = Gam9thet;
% spring
GamThet(M4+1:M5,M4+1:M5) = Gam5thetu;
GamThet(M4+1:M5,M9+M4+1:M9+M5) = Gam5thetv;
GamThet(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = Gam5thetw;
GamThet(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = Gam5thetthet;

Gam = [GamU; GamV; GamW; GamThet];
% matrix corresponding to: [u1 -> u9; v1 -> v9; w1 -> w9; thet1 -> thet9]'


%% DAMPING MATRICES %%
R1u = psi(1)*Gam1u - (sigA * eye(size(Gam1u)));
R2u = psi(2)*Gam2u - (sigA * eye(size(Gam2u)));
R3u = psi(3)*Gam3u - (sigA * eye(size(Gam3u)));
R4u = psi(4)*Gam4u - (sigA * eye(size(Gam4u)));

R5uu = psi(5)*Gam5uu - (sigA * eye(size(Gam5uu)));
R5uv = psi(5)*Gam5uv;
R5uw = psi(5)*Gam5uw;
R5uthet = psi(5)*Gam5uthet;

R6u = psi(6)*Gam6u - (sigA * eye(size(Gam6u)));
R7u = psi(7)*Gam7u - (sigA * eye(size(Gam7u)));
R8u = psi(8)*Gam8u - (sigA * eye(size(Gam8u)));
R9u = psi(9)*Gam9u - (sigA * eye(size(Gam9u)));


R1v = psi(1)*Gam1v - (sigA * eye(size(Gam1v)));
R2v = psi(2)*Gam2v - (sigA * eye(size(Gam2v)));
R3v = psi(3)*Gam3v - (sigA * eye(size(Gam3v)));
R4v = psi(4)*Gam4v - (sigA * eye(size(Gam4v)));

R5vu = psi(5)*Gam5vu;
R5vv = psi(5)*Gam5vv - (sigA * eye(size(Gam5vv)));
R5vw = psi(5)*Gam5vw;
R5vthet = psi(5)*Gam5vthet;

R6v = psi(6)*Gam6v - (sigA * eye(size(Gam6v)));
R7v = psi(7)*Gam7v - (sigA * eye(size(Gam7v)));
R8v = psi(8)*Gam8v - (sigA * eye(size(Gam8v)));
R9v = psi(9)*Gam9v - (sigA * eye(size(Gam9v)));

R1w = psi(1)*Gam1w - (sigB * eye(size(Gam1w)));
R2w = psi(2)*Gam2w - (sigB * eye(size(Gam2w)));
R3w = psi(3)*Gam3w - (sigB * eye(size(Gam3w)));
R4w = psi(4)*Gam4w - (sigB * eye(size(Gam4w)));

R5wu = psi(5)*Gam5wu;
R5wv = psi(5)*Gam5wv;
R5ww = psi(5)*Gam5ww - (sigB * eye(size(Gam5ww)));
R5wthet = psi(5)*Gam5wthet;

R6w = psi(6)*Gam6w - (sigB * eye(size(Gam6w)));
R7w = psi(7)*Gam7w - (sigB * eye(size(Gam7w)));
R8w = psi(8)*Gam8w - (sigB * eye(size(Gam8w)));
R9w = psi(9)*Gam9w - (sigB * eye(size(Gam9w)));


R1thet = psi(1)*Gam1thet - (sigC * eye(size(Gam1thet)));
R2thet = psi(2)*Gam2thet - (sigC * eye(size(Gam2thet)));
R3thet = psi(3)*Gam3thet - (sigC * eye(size(Gam3thet)));
R4thet = psi(4)*Gam4thet - (sigC * eye(size(Gam4thet)));

R5thetu = psi(5)*Gam5thetu;
R5thetv = psi(5)*Gam5thetv;
R5thetw = psi(5)*Gam5thetw;
R5thetthet = psi(5)*Gam5thetthet - (sigC * eye(size(Gam5thetthet)));

R6thet = psi(6)*Gam6thet - (sigC * eye(size(Gam6thet)));
R7thet = psi(7)*Gam7thet - (sigC * eye(size(Gam7thet)));
R8thet = psi(8)*Gam8thet - (sigC * eye(size(Gam8thet)));
R9thet = psi(9)*Gam9thet - (sigC * eye(size(Gam9thet)));


%% CONSOLIDATE DAMPING MATRICES %%
% u polarisation for all segments
RU = zeros(M9,4*M9);
RU(1:M1,1:M1) = R1u;
RU(M1+1:M2,M1+1:M2) = R2u;
RU(M2+1:M3,M2+1:M3) = R3u;
RU(M3+1:M4,M3+1:M4) = R4u;

RU(M4+1:M5,M4+1:M5) = R5uu;
RU(M4+1:M5,M9+M4+1:M9+M5) = R5uv;
RU(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = R5uw;
RU(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = R5uthet;

RU(M5+1:M6,M5+1:M6) = R6u;
RU(M6+1:M7,M6+1:M7) = R7u;
RU(M7+1:M8,M7+1:M8) = R8u;
RU(M8+1:M9,M8+1:M9) = R9u;

% v polarisation for all segments
RV = zeros(M9,4*M9);
RV(1:M1,M9+(1:M1)) = R1v;
RV(M1+1:M2,M9+(M1+1:M2)) = R2v;
RV(M2+1:M3,M9+(M2+1:M3)) = R3v;
RV(M3+1:M4,M9+(M3+1:M4)) = R4v;

RV(M4+1:M5,M4+1:M5) = R5vu;
RV(M4+1:M5,M9+M4+1:M9+M5) = R5vv;
RV(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = R5vw;
RV(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = R5vthet;

RV(M5+1:M6,M9+(M5+1:M6)) = R6v;
RV(M6+1:M7,M9+(M6+1:M7)) = R7v;
RV(M7+1:M8,M9+(M7+1:M8)) = R8v;
RV(M8+1:M9,M9+(M8+1:M9)) = R9v;

% w polarisation for all segments
RW = zeros(M9,4*M9);
RW(1:M1,(2*M9)+(1:M1)) = R1w;
RW(M1+1:M2,(2*M9)+(M1+1:M2)) = R2w;
RW(M2+1:M3,(2*M9)+(M2+1:M3)) = R3w;
RW(M3+1:M4,(2*M9)+(M3+1:M4)) = R4w;

RW(M4+1:M5,M4+1:M5) = R5wu;
RW(M4+1:M5,M9+M4+1:M9+M5) = R5wv;
RW(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = R5ww;
RW(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = R5wthet;

RW(M5+1:M6,(2*M9)+(M5+1:M6)) = R6w;
RW(M6+1:M7,(2*M9)+(M6+1:M7)) = R7w;
RW(M7+1:M8,(2*M9)+(M7+1:M8)) = R8w;
RW(M8+1:M9,(2*M9)+(M8+1:M9)) = R9w;

% thet polarisation for all segments
RThet = zeros(M9,4*M9);
RThet(1:M1,(3*M9)+(1:M1)) = R1thet;
RThet(M1+1:M2,(3*M9)+(M1+1:M2)) = R2thet;
RThet(M2+1:M3,(3*M9)+(M2+1:M3)) = R3thet;
RThet(M3+1:M4,(3*M9)+(M3+1:M4)) = R4thet;

RThet(M4+1:M5,M4+1:M5) = R5thetu;
RThet(M4+1:M5,M9+M4+1:M9+M5) = R5thetv;
RThet(M4+1:M5,(2*M9)+M4+1:(2*M9)+M5) = R5thetw;
RThet(M4+1:M5,(3*M9)+M4+1:(3*M9)+M5) = R5thetthet;

RThet(M5+1:M6,(3*M9)+(M5+1:M6)) = R6thet;
RThet(M6+1:M7,(3*M9)+(M6+1:M7)) = R7thet;
RThet(M7+1:M8,(3*M9)+(M7+1:M8)) = R8thet;
RThet(M8+1:M9,(3*M9)+(M8+1:M9)) = R9thet;

R = [RU; RV; RW; RThet];
% matrix corresponding to: [u1 -> u9; v1 -> v9; w1 -> w9; thet1 -> thet9]'


%% CONVENIENT INDEXING %%
% boundary conditions
u1LB = [1 2]; u9RB = [M9-1 M9];

v1LB = M9 + u1LB; v9RB = M9 + u9RB;
w1LB = (2*M9) + u1LB(1); w9RB = (2*M9) + u9RB(2);
thet1LB = (3*M9) + u1LB(1); thet9RB = (3*M9) + u9RB(2);

% consolidate boundary indexing
uB = [u1LB u9RB]; vB = [v1LB v9RB];
wB = [w1LB w9RB]; thetB = [thet1LB thet9RB];

% continuity conditions
uC12 = [M1-1 M1 M1+1 M1+2];
uC23 = [M2-1 M2 M2+1 M2+2];
uC34 = [M3-1 M3 M3+1 M3+2];
uC45 = [M4-1 M4 M4+1 M4+2];
uC56 = [M5-1 M5 M5+1 M5+2];
uC67 = [M6-1 M6 M6+1 M6+2];
uC78 = [M7-1 M7 M7+1 M7+2];
uC89 = [M8-1 M8 M8+1 M8+2];

vC12 = M9 + uC12;
vC23 = M9 + uC23;
vC34 = M9 + uC34;
vC45 = M9 + uC45;
vC56 = M9 + uC56;
vC67 = M9 + uC67;
vC78 = M9 + uC78;
vC89 = M9 + uC89;

wC12 = (2*M9) + uC12(2:3);
wC23 = (2*M9) + uC23(2:3);
wC34 = (2*M9) + uC34(2:3);
wC45 = (2*M9) + uC45(2:3);
wC56 = (2*M9) + uC56(2:3);
wC67 = (2*M9) + uC67(2:3);
wC78 = (2*M9) + uC78(2:3);
wC89 = (2*M9) + uC89(2:3);

thetC12 = (3*M9) + uC12(2:3);
thetC23 = (3*M9) + uC23(2:3);
thetC34 = (3*M9) + uC34(2:3);
thetC45 = (3*M9) + uC45(2:3);
thetC56 = (3*M9) + uC56(2:3);
thetC67 = (3*M9) + uC67(2:3);
thetC78 = (3*M9) + uC78(2:3);
thetC89 = (3*M9) + uC89(2:3);

% consolidate continuity indexing
uC = [uC12 uC23 uC34 uC45 uC56 uC67 uC78 uC89];
vC = [vC12 vC23 vC34 vC45 vC56 vC67 vC78 vC89];
wC = [wC12 wC23 wC34 wC45 wC56 wC67 wC78 wC89];
thetC = [thetC12 thetC23 thetC34 thetC45 thetC56 thetC67 thetC78 thetC89];

% start and end points
u1L = 1; u1R = M1;
v1L = u1L + M9; v1R = u1R + M9;
w1L = u1L + (2*M9); w1R = u1R + (2*M9);
thet1L = u1L + (3*M9); thet1R = u1R + (3*M9);

u2L = M1+1; u2R = M2;
v2L = u2L + M9; v2R = u2R + M9;
w2L = u2L + (2*M9); w2R = u2R + (2*M9);
thet2L = u2L + (3*M9); thet2R = u2R + (3*M9);

u3L = M2+1; u3R = M3;
v3L = u3L + M9; v3R = u3R + M9;
w3L = u3L + (2*M9); w3R = u3R + (2*M9);
thet3L = u3L + (3*M9); thet3R = u3R + (3*M9);

u4L = M3+1; u4R = M4;
v4L = u4L + M9; v4R = u4R + M9;
w4L = u4L + (2*M9); w4R = u4R + (2*M9);
thet4L = u4L + (3*M9); thet4R = u4R + (3*M9);

u5L = M4+1; u5R = M5;
v5L = u5L + M9; v5R = u5R + M9;
w5L = u5L + (2*M9); w5R = u5R + (2*M9);
thet5L = u5L + (3*M9); thet5R = u5R + (3*M9);

u6L = M5+1; u6R = M6;
v6L = u6L + M9; v6R = u6R + M9;
w6L = u6L + (2*M9); w6R = u6R + (2*M9);
thet6L = u6L + (3*M9); thet6R = u6R + (3*M9);

u7L = M6+1; u7R = M7;
v7L = u7L + M9; v7R = u7R + M9;
w7L = u7L + (2*M9); w7R = u7R + (2*M9);
thet7L = u7L + (3*M9); thet7R = u7R + (3*M9);

u8L = M7+1; u8R = M8;
v8L = u8L + M9; v8R = u8R + M9;
w8L = u8L + (2*M9); w8R = u8R + (2*M9);
thet8L = u8L + (3*M9); thet8R = u8R + (3*M9);

u9L = M8+1; u9R = M9;
v9L = u9L + M9; v9R = u9R + M9;
w9L = u9L + (2*M9); w9R = u9R + (2*M9);
thet9L = u9L + (3*M9); thet9R = u9R + (3*M9);

% interior domains
u1i = (u1L+2):(u1R-2); u2i = (u2L+2):(u2R-2); u3i = (u3L+2):(u3R-2);
u4i = (u4L+2):(u4R-2); u5i = (u5L+2):(u5R-2); u6i = (u6L+2):(u6R-2);
u7i = (u7L+2):(u7R-2); u8i = (u8L+2):(u8R-2); u9i = (u9L+2):(u9R-2);

v1i = (v1L+2):(v1R-2); v2i = (v2L+2):(v2R-2); v3i = (v3L+2):(v3R-2);
v4i = (v4L+2):(v4R-2); v5i = (v5L+2):(v5R-2); v6i = (v6L+2):(v6R-2);
v7i = (v7L+2):(v7R-2); v8i = (v8L+2):(v8R-2); v9i = (v9L+2):(v9R-2);

w1i = (w1L+1):(w1R-1); w2i = (w2L+1):(w2R-1); w3i = (w3L+1):(w3R-1);
w4i = (w4L+1):(w4R-1); w5i = (w5L+1):(w5R-1); w6i = (w6L+1):(w6R-1);
w7i = (w7L+1):(w7R-1); w8i = (w8L+1):(w8R-1); w9i = (w9L+1):(w9R-1);

thet1i = (thet1L+1):(thet1R-1); thet2i = (thet2L+1):(thet2R-1); thet3i = (thet3L+1):(thet3R-1);
thet4i = (thet4L+1):(thet4R-1); thet5i = (thet5L+1):(thet5R-1); thet6i = (thet6L+1):(thet6R-1);
thet7i = (thet7L+1):(thet7R-1); thet8i = (thet8L+1):(thet8R-1); thet9i = (thet9L+1):(thet9R-1);

% consolidate interior domains
ui = [u1i u2i u3i u4i u5i u6i u7i u8i u9i];
vi = [v1i v2i v3i v4i v5i v6i v7i v8i v9i];
wi = [w1i w2i w3i w4i w5i w6i w7i w8i w9i];
theti = [thet1i thet2i thet3i thet4i thet5i thet6i thet7i thet8i thet9i];


%% ENFORCE BOUNDARY CONDITIONS %%
% zero out rows to be replaced
Gam(u1LB,u1L:u1R) = 0; Gam(u9RB,u9L:u9R) = 0;
Gam(v1LB,v1L:v1R) = 0; Gam(v9RB,v9L:v9R) = 0;
Gam(w1LB,w1L:w1R) = 0; Gam(w9RB,w9L:w9R) = 0;
Gam(thet1LB,thet1L:thet1R) = 0; Gam(thet9RB,thet9L:thet9R) = 0;

% zero them for the damping matrix too
R(u1LB,u1L:u1R) = 0; R(u9RB,u9L:u9R) = 0;
R(v1LB,v1L:v1R) = 0; R(v9RB,v9L:v9R) = 0;
R(w1LB,w1L:w1R) = 0; R(w9RB,w9L:w9R) = 0;
R(thet1LB,thet1L:thet1R) = 0; R(thet9RB,thet9L:thet9R) = 0;

if f_BCs == "pinned"
    % LBCs
    Gam(u1LB(1),u1LB(1)) = 1;
    Gam(u1LB(2),u1L:u1R) = D12(1,:);

    Gam(v1LB(1),v1LB(1)) = 1;
    Gam(v1LB(2),v1L:v1R) = D12(1,:);

    Gam(w1LB,w1LB) = 1;

    Gam(thet1LB,thet1LB) = 1;

    % RBCs
    Gam(u9RB(2),u9RB(2)) = 1;
    Gam(u9RB(1),u9L:u9R) = D92(N9,:);

    Gam(v9RB(2),v9RB(2)) = 1;
    Gam(v9RB(1),v9L:v9R) = D92(N9,:);

    Gam(w9RB,w9RB) = 1;

    Gam(thet9RB,thet9RB) = 1;

elseif f_BCs == "clamped"
    % LBCs
    Gam(u1LB(1),u1LB(1)) = 1;
    Gam(u1LB(2),u1L:u1R) = D11(1,:);

    Gam(v1LB(1),v1LB(1)) = 1;
    Gam(v1LB(2),v1L:v1R) = D11(1,:);

    Gam(w1LB,w1LB) = 1;

    Gam(thet1LB,thet1LB) = 1;

    % RBCs
    Gam(u9RB(2),u9RB(2)) = 1;
    Gam(u9RB(1),u9L:u9R) = D91(N9,:);

    Gam(v9RB(2),v9RB(2)) = 1;
    Gam(v9RB(1),v9L:v9R) = D91(N9,:);

    Gam(w9RB,w9RB) = 1;

    Gam(thet9RB,thet9RB) = 1;
end


%% CONTINUITY: SYSTEMS 1 -> 2 %%
% zero out rows to be replaced
Gam(uC12,u1L:u2R) = 0;
Gam(vC12,v1L:v2R) = 0;
Gam(wC12,w1L:w2R) = 0;
Gam(thetC12,thet1L:thet2R) = 0;

% zero for damping matrix too
R(uC12,u1L:u2R) = 0;
R(vC12,v1L:v2R) = 0;
R(wC12,w1L:w2R) = 0;
R(thetC12,thet1L:thet2R) = 0;

% u
Gam(uC12(1),u1R) = 1;
Gam(uC12(1),u2L) = -1;

Gam(uC12(2),u1L:u1R) = D11(N1,:);
Gam(uC12(2),u2L:u2R) = -D21(1,:);

Gam(uC12(3),u1L:u1R) = E(1)*I(1) * D12(N1,:);
Gam(uC12(3),u2L:u2R) = -E(2)*I(2) * D22(1,:);

Gam(uC12(4),u1L:u1R) = E(1)*I(1) * D13(N1,:);
Gam(uC12(4),u2L:u2R) = -E(2)*I(2) * D23(1,:);

% v
Gam(vC12(1),v1R) = 1;
Gam(vC12(1),v2L) = -1;

Gam(vC12(2),v1L:v1R) = D11(N1,:);
Gam(vC12(2),v2L:v2R) = -D21(1,:);

Gam(vC12(3),v1L:v1R) = E(1)*I(1) * D12(N1,:);
Gam(vC12(3),v2L:v2R) = -E(2)*I(2) * D22(1,:);

Gam(vC12(4),v1L:v1R) = E(1)*I(1) * D13(N1,:);
Gam(vC12(4),v2L:v2R) = -E(2)*I(2) * D23(1,:);

% w
Gam(wC12(1),w1R) = 1;
Gam(wC12(1),w2L) = -1;

Gam(wC12(2),w1L:w1R) = E(1)*CSA(1) * D11(N1,:);
Gam(wC12(2),w2L:w2R) = -E(2)*CSA(2) * D21(1,:);

% theta_w
Gam(thetC12(1),thet1R) = 1;
Gam(thetC12(1),thet2L) = -1;

Gam(thetC12(2),thet1L:thet1R) = 2*G(1)*I(1) * D11(N1,:);
Gam(thetC12(2),thet2L:thet2R) = -2*G(2)*I(2) * D21(1,:);


%% CONTINUITY: SYSTEMS 2 -> 3 %%
% zero out rows to be replaced
Gam(uC23,u2L:u3R) = 0;
Gam(vC23,v2L:v3R) = 0;
Gam(wC23,w2L:w3R) = 0;
Gam(thetC23,thet2L:thet3R) = 0;

% zero for damping matrix too
R(uC23,u2L:u3R) = 0;
R(vC23,v2L:v3R) = 0;
R(wC23,w2L:w3R) = 0;
R(thetC23,thet2L:thet3R) = 0;

% u
Gam(uC23(1),u2R) = 1;
Gam(uC23(1),u3L) = -1;

Gam(uC23(2),u2L:u2R) = D21(N2,:);
Gam(uC23(2),u3L:u3R) = -D31(1,:);

Gam(uC23(3),u2L:u2R) = E(2)*I(2) * D22(N2,:);
Gam(uC23(3),u3L:u3R) = -E(3)*I(3) * D32(1,:);

Gam(uC23(4),u2L:u2R) = E(2)*I(2) * D23(N2,:);
Gam(uC23(4),u3L:u3R) = -E(3)*I(3) * D33(1,:);

% v
Gam(vC23(1),v2R) = 1;
Gam(vC23(1),v3L) = -1;

Gam(vC23(2),v2L:v2R) = D21(N2,:);
Gam(vC23(2),v3L:v3R) = -D31(1,:);

Gam(vC23(3),v2L:v2R) = E(2)*I(2) * D22(N2,:);
Gam(vC23(3),v3L:v3R) = E(3)*I(3) * -D32(1,:);

Gam(vC23(4),v2L:v2R) = E(2)*I(2) * D23(N2,:);
Gam(vC23(4),v3L:v3R) = E(3)*I(3) * -D33(1,:);

% w
Gam(wC23(1),w2R) = 1;
Gam(wC23(1),w3L) = -1;

Gam(wC23(2),w2L:w2R) = E(2)*CSA(2) * D21(N2,:);
Gam(wC23(2),w3L:w3R) = E(3)*CSA(3) * -D31(1,:);

% theta_w
Gam(thetC23(1),thet2R) = 1;
Gam(thetC23(1),thet3L) = -1;

Gam(thetC23(2),thet2L:thet2R) = 2*G(2)*I(2) * D21(N2,:);
Gam(thetC23(2),thet3L:thet3R) = 2*G(3)*I(3) * -D31(1,:);


%% CONTINUITY: SYSTEMS 3 -> 4 %%
% zero out rows to be replaced
Gam(uC34,:) = 0;
Gam(vC34,:) = 0;
Gam(wC34,:) = 0;
Gam(thetC34,:) = 0;

% zero for damping matrix too
R(uC34,:) = 0;
R(vC34,:) = 0;
R(wC34,:) = 0;
R(thetC34,:) = 0;

% u3 - u4 = 0
Gam(uC34(1),u3R) = 1;
Gam(uC34(1),u4L) = -1;
% thetV3 + thetW4 = 0
Gam(uC34(2),u3L:u3R) = D31(N3,:);
Gam(uC34(2),thet4L) = 1; % thetW4
% mV3 + mW4 = 0
Gam(uC34(3),u3L:u3R) = E(3)*I(3) * D32(N3,:); % mV3
Gam(uC34(3),thet4L:thet4R) = 2*G(4)*I(4) * D41(1,:); % mW4
% pU3 - pU4 = 0
Gam(uC34(4),u3L:u3R) = -E(3)*I(3) * D33(N3,:); % pU3
Gam(uC34(4),u4L:u4R) = E(4)*I(4) * D43(1,:); % pU4


% v3 + w4 = 0
Gam(vC34(1),v3R) = 1;
Gam(vC34(1),w4L) = 1; 
% thetU3 - thetU4 = 0
Gam(vC34(2),v3L:v3R) = -D31(N3,:); % thetU3
Gam(vC34(2),v4L:v4R) = D41(1,:); % thetU4
% mU3 - mU4 = 0
Gam(vC34(3),v3L:v3R) = -E(3)*I(3) * D32(N3,:); % mU3
Gam(vC34(3),v4L:v4R) = E(4)*I(4) * D42(1,:); % mU4
% pV3 + pW4 = 0
Gam(vC34(4),v3L:v3R) = -E(3)*I(3) * D33(N3,:); % pV3
Gam(vC34(4),w4L:w4R) = E(4)*CSA(4) * D41(1,:); % pW4

% w3 - v4 = 0
Gam(wC34(1),w3R) = 1;
Gam(wC34(1),v4L) = -1;
% pW3 - pV4 = 0
Gam(wC34(2),w3L:w3R) = E(3)*CSA(3) * D31(N3,:); % pW3
Gam(wC34(2),v4L:v4R) = E(4)*I(4) * D43(1,:); % pV4

% thetW3 - thetV4 = 0
Gam(thetC34(1),thet3R) = 1; % thetW3
Gam(thetC34(1),u4L:u4R) = -D41(1,:); % thetV4
% mW3 - mV4 = 0
Gam(thetC34(2),thet3L:thet3R) = 2*G(3)*I(3) * D31(N3,:); % mW3
Gam(thetC34(2),u4L:u4R) = -E(4)*I(4) * D42(1,:); % mV4


%% CONTINUITY: SYSTEMS 4 -> 5  %%
% zero out rows to be replaced
Gam(uC45,:) = 0;
Gam(vC45,:) = 0;
Gam(wC45,:) = 0;
Gam(thetC45,:) = 0;

% zero for damping matrix too
R(uC45,:) = 0;
R(vC45,:) = 0;
R(wC45,:) = 0;
R(thetC45,:) = 0;

% REMEMBER: A = sin(alp), B = cos(alp)

% u4R - B*w5L + A*v5L = 0
Gam(uC45(1),u4R) = 1; % u4
Gam(uC45(1),w5L) = -B; % w5
Gam(uC45(1),v5L) = A; % v5

% thetV4R - B*thetV5L - A*thetW5L = 0
Gam(uC45(2),u4L:u4R) = D41(N4,:); % thetV4

Gam(uC45(2),u5L:u5R) = -B*D51(1,:); % thetV5 u
Gam(uC45(2),v5L) = B*kap5*mu5; % thetV5 v
Gam(uC45(2),w5L) = -B*kap5; % thetV5 w

Gam(uC45(2),thet5L) = -A; % thetW5

% mV4R - B*mV5L - A*mW5L = 0
Gam(uC45(3),u4L:u4R) = E(4)*I(4) * D42(N4,:); % mV4

Gam(uC45(3),u5L:u5R) = -B*E(5)*I(5) * D52(1,:); % mV5: u
Gam(uC45(3),u5L) = Gam(uC45(3),u5L) + (B*E(5)*I(5) * (kap5^2)*(mu5^2)); % mV5: u
Gam(uC45(3),v5L:v5R) = 2*B*E(5)*I(5)*kap5*mu5 * D51(1,:); % mV5: v
Gam(uC45(3),w5L:w5R) = -B*E(5)*I(5)*kap5 * D51(1,:); % mV5: w

Gam(uC45(3),u5L) = Gam(uC45(3),u5L) + (-2*A*G(5)*I(5) * (kap5^2)*mu5); % mW5: u
Gam(uC45(3),v5L:v5R) = Gam(uC45(3),v5L:v5R) + (-2*A*G(5)*I(5) * kap5 * D51(1,:)); % mW5: v
Gam(uC45(3),thet5L:thet5R) = -2*A*G(5)*I(5) * D51(1,:); % mW5: thetW


% pU4R - B*pW5L + A*pV5L = 0
Gam(uC45(4),u4L:u4R) = -E(4)*I(4) * D43(N4,:); % pU4

Gam(uC45(4),u5L) = B*E(5)*CSA(5) * kap5; % pW5 u
Gam(uC45(4),w5L:w5R) = -B*E(5)*CSA(5) * D51(1,:); % pW5 w

Gam(uC45(4),u5L:u5R) = Gam(uC45(4),u5L:u5R) + (-3*A*E(5)*I(5) * kap5*mu5 * D52(1,:)); % pV5 u
Gam(uC45(4),u5L) = Gam(uC45(4),u5L) + (A*E(5)*I(5) * (kap5^3)*(mu5^3)) + (2*A*G(5)*I(5) * (kap5^3)*mu5); % pV5 u
Gam(uC45(4),v5L:v5R) = (A*E(5)*I(5) * ((3*(kap5^2)*(mu5^2)*D51(1,:)) - D53(1,:))) + (2*A*G(5)*I(5) * (kap5^2)*D51(1,:)); % pV5 v
Gam(uC45(4),w5L:w5R) = Gam(uC45(4),w5L:w5R) + (-2*A*E(5)*I(5) * (kap5^2)*mu5 * D51(1,:)); % pV5 w
Gam(uC45(4),thet5L:thet5R) = (A*E(5)*I(5) * kap5*D51(1,:)) + (2*A*G(5)*I(5) * kap5*D51(1,:)); % pV5 thetW

% v4R - B*v5L - A*w5L = 0
Gam(vC45(1),v4R) = 1; % v4
Gam(vC45(1),v5L) = -B; % v5
Gam(vC45(1),w5L) = -A; % w5

% thetU4R - B*thetW5L + A*thetV5L = 0
Gam(vC45(2),v4L:v4R) = -D41(N4,:); % thetU4

Gam(vC45(2),thet5L) = -B; % thetW5

Gam(vC45(2),u5L:u5R) = A*D51(1,:); % thetV5 u
Gam(vC45(2),v5L) = -A*kap5*mu5; % thetV5 v
Gam(vC45(2),w5L) = A*kap5; % thetV5 w

% mU4R - B*mW5L + A*mV5L = 0
Gam(vC45(3),v4L:v4R) = -E(4)*I(4) * D42(N4,:); % mU4

Gam(vC45(3),u5L) = -B*2*G(5)*I(5) * (kap5^2)*mu5; % mW5 u 
Gam(vC45(3),v5L:v5R) = -B*2*G(5)*I(5) * kap5 * D51(1,:); % mW5 v
Gam(vC45(3),thet5L:thet5R) = -B*2*G(5)*I(5) * D51(1,:); % mW5 thetW

Gam(vC45(3),u5L:u5R) = Gam(vC45(3),u5L:u5R) + (A*E(5)*I(5) * D52(1,:)); % mV5: u
Gam(vC45(3),u5L) = Gam(vC45(3),u5L) - (A*E(5)*I(5) * (kap5^2)*(mu5^2)); % mV5: u
Gam(vC45(3),v5L:v5R) = Gam(vC45(3),v5L:v5R) + (-2*A*E(5)*I(5)*kap5*mu5 * D51(1,:)); % mV5: v
Gam(vC45(3),w5L:w5R) = A*E(5)*I(5)*kap5 * D51(1,:); % mV5: w


% pV4R - B*pV5L - A*pW5L = 0
Gam(vC45(4),v4L:v4R) = -E(4)*I(4) * D43(N4,:); % pV4

Gam(vC45(4),u5L:u5R) = 3*B*E(5)*I(5) * kap5*mu5 * D52(1,:); % pV5 u
Gam(vC45(4),u5L) = Gam(vC45(4),u5L) - (B*E(5)*I(5) * (kap5^3)*(mu5^3)) - (2*B*G(5)*I(5) * (kap5^3)*mu5); % pV5 u
Gam(vC45(4),v5L:v5R) = -(B*E(5)*I(5) * ((3*(kap5^2)*(mu5^2)*D51(1,:)) - D53(1,:))) - (2*B*G(5)*I(5) * (kap5^2)*D51(1,:)); % pV5 v
Gam(vC45(4),w5L:w5R) = 2*B*E(5)*I(5) * (kap5^2)*mu5 * D51(1,:); % pV5 w
Gam(vC45(4),thet5L:thet5R) = -(B*E(5)*I(5) * kap5*D51(1,:)) - (2*B*G(5)*I(5) * kap5*D51(1,:)); % pV5 thetW

Gam(vC45(4),u5L) =  Gam(vC45(4),u5L) + (A*E(5)*CSA(5) * kap5); % pW5 u
Gam(vC45(4),w5L:w5R) = Gam(vC45(4),w5L:w5R) + (-A*E(5)*CSA(5) * D51(1,:)); % pW5 w

% w4R + u5L = 0
Gam(wC45(1),w4R) = 1;
Gam(wC45(1),u5L) = 1;

% pW4R + pU5L = 0
Gam(wC45(2),w4L:w4R) = E(4)*CSA(4) * D41(N4,:); % pW4

Gam(wC45(2),u5L:u5R) = E(5)*I(5) * ((3*(kap5^2)*(mu5^2)*D51(1,:)) - D53(1,:)); % pU5 u
Gam(wC45(2),v5L:v5R) = E(5)*I(5) * (3 * kap5 * mu5 * D52(1,:)); % pU5 v
Gam(wC45(2),v5L) = Gam(wC45(2),v5L) - (E(5)*I(5) * ((kap5^3)*(mu5^3))); % pU5 v
Gam(wC45(2),w5L:w5R) = -E(5)*I(5) * (kap5*D52(1,:)); % pU5 w
Gam(wC45(2),w5L) = Gam(wC45(2),w5L) + (E(5)*I(5) * ((kap5^3)*(mu5^2))); % pU5 w
Gam(wC45(2),thet5L) = -E(5)*I(5) * ((kap5^2)*mu5); % pU5 thetW

% thetW4R + thetU5L = 0
Gam(thetC45(1),thet4R) = 1; % thetW4

Gam(thetC45(1),u5L) = -kap5*mu5; % thetU5 u
Gam(thetC45(1),v5L:v5R) = -D51(1,:); % thetU5 v

% mW4R + mU5L = 0
Gam(thetC45(2),thet4L:thet4R) = 2*G(4)*I(4) * D41(N4,:); % mW4

Gam(thetC45(2),u5L:u5R) = -2*E(5)*I(5) * (kap5*mu5) * D51(1,:); % mU4 u
Gam(thetC45(2),v5L:v5R) = -E(5)*I(5) * D52(1,:); % mU4 v
Gam(thetC45(2),v5L) = Gam(thetC45(2),v5L) + (E(5)*I(5) * ((kap5^2)*(mu5^2))); % mU4 v
Gam(thetC45(2),w5L) = -E(5)*I(5) * ((kap5^2)*mu5); % mU4 w
Gam(thetC45(2),thet5L) = E(5)*I(5) * kap5; % mU4 thetW


%% CONTINUITY: SYSTEMS 5 -> 6 %%
% zero out rows to be replaced
Gam(uC56,:) = 0;
Gam(vC56,:) = 0;
Gam(wC56,:) = 0;
Gam(thetC56,:) = 0;

% zero for damping matrix too
R(uC56,:) = 0;
R(vC56,:) = 0;
R(wC56,:) = 0;
R(thetC56,:) = 0;

% u5 - w6 = 0
Gam(uC56(1),w6L) = -1; % w6
Gam(uC56(1),u5R) = 1; % u5

% -B*thetV5 - A*thetW5 + thetV6 = 0
Gam(uC56(2),u6L:u6R) = D61(1,:); % thetV6

Gam(uC56(2),u5L:u5R) = -B*D51(N5,:); % thetV5 u
Gam(uC56(2),v5R) = B*kap5*mu5; % thetV5 v
Gam(uC56(2),w5R) = -B*kap5; % thetV5 w

Gam(uC56(2),thet5R) = -A; % thetW5

% -B*mV5 - A*mW5 + mV6 = 0
Gam(uC56(3),u6L:u6R) = E(6)*I(6) * D62(1,:); % mV6

Gam(uC56(3),u5L:u5R) = -B*E(5)*I(5) * D52(N5,:); % mV5: u
Gam(uC56(3),u5R) = Gam(uC56(3),u5R) + (B*E(5)*I(5) * (kap5^2)*(mu5^2)); % mV5: u
Gam(uC56(3),v5L:v5R) = 2*B*E(5)*I(5)*kap5*mu5 * D51(N5,:); % mV5: v
Gam(uC56(3),w5L:w5R) = -B*E(5)*I(5)*kap5 * D51(N5,:); % mV5: w

Gam(uC56(3),u5R) = Gam(uC56(3),u5R) - (2*A*G(5)*I(5) * (kap5^2)*mu5); % mW5: u
Gam(uC56(3),v5L:v5R) = Gam(uC56(3),v5L:v5R) - (2*A*G(5)*I(5) * kap5 * D51(N5,:)); % mW5: v
Gam(uC56(3),thet5L:thet5R) = -2*A*G(5)*I(5) * D51(N5,:); % mW5: thetW

% pU5 - pW6 = 0
Gam(uC56(4),w6L:w6R) = -E(6)*CSA(6) * D61(1,:); % pW6

Gam(uC56(4),u5L:u5R) = E(5)*I(5) * ((3*(kap5^2)*(mu5^2)*D51(N5,:)) - D53(N5,:)); % pU5 u
Gam(uC56(4),v5L:v5R) = E(5)*I(5) *  (3 * kap5 * mu5 * D52(N5,:)); % pU5 v
Gam(uC56(4),v5R) = Gam(uC56(4),v5R) - (E(5)*I(5) * ((kap5^3)*(mu5^3))); % pU5 v
Gam(uC56(4),w5L:w5R) = -E(5)*I(5)*kap5 * D52(N5,:); % pU5 w
Gam(uC56(4),w5R) = Gam(uC56(4),w5R) + (E(5)*I(5) * ((kap5^3)*(mu5^2))); % pU5 w
Gam(uC56(4),thet5R) = -E(5)*I(5) * ((kap5^2)*mu5); % pU5 thetW

% B*w5 - A*v5 + u6 = 0
Gam(vC56(1),u6L) = 1; % u6

Gam(vC56(1),w5R) = B; % w5
Gam(vC56(1),v5R) = -A; % v5

% thetU5 - thetW6 = 0
Gam(vC56(2),thet6L) = -1; % thetW6

Gam(vC56(2),u5R) = -kap5*mu5; % thetU5 u
Gam(vC56(2),v5L:v5R) = -D51(N5,:); % thetU5 v

% mU5 - mW6 = 0
Gam(vC56(3),thet6L:thet6R) = -2*G(6)*I(6) * D61(1,:); % mW6

Gam(vC56(3),u5L:u5R) = -2*E(5)*I(5) * (kap5*mu5) * D51(N5,:); % mU5 u
Gam(vC56(3),v5L:v5R) = -E(5)*I(5) * D52(N5,:); % mU5 v
Gam(vC56(3),v5R) = Gam(vC56(3),v5R) + (E(5)*I(5) * ((kap5^2)*(mu5^2))); % mU5 v
Gam(vC56(3),w5R) = -E(5)*I(5) * ((kap5^2)*mu5); % mU5 w
Gam(vC56(3),thet5R) = E(5)*I(5) * kap5; % mU5 thetW

% -B*pV5 - A*pW5 + pV6 = 0
Gam(vC56(4),v6L:v6R) = -E(6)*I(6) * D63(1,:); % pV6

Gam(vC56(4),u5L:u5R) = 3*B*E(5)*I(5) * kap5*mu5 * D52(N5,:); % pV5 u
Gam(vC56(4),u5R) = Gam(vC56(4),u5R) - (B*E(5)*I(5) * (kap5^3)*(mu5^3)) - (2*B*G(5)*I(5) * (kap5^3)*mu5); % pV5 u
Gam(vC56(4),v5L:v5R) = -(B*E(5)*I(5) * ((3*(kap5^2)*(mu5^2)*D51(N5,:)) - D53(N5,:))) - (2*B*G(5)*I(5) * (kap5^2)*D51(N5,:)); % pV5 v
Gam(vC56(4),w5L:w5R) = 2*B*E(5)*I(5) * (kap5^2)*mu5 * D51(N5,:); % pV5 w
Gam(vC56(4),thet5L:thet5R) = -(B*E(5)*I(5) * kap5*D51(N5,:)) - (2*B*G(5)*I(5) * kap5*D51(N5,:)); % pV5 thetW

Gam(vC56(4),u5R) = Gam(vC56(4),u5R) + (A*E(5)*CSA(5) * kap5); % pW5 u
Gam(vC56(4),w5L:w5R) = Gam(vC56(4),w5L:w5R) - (A*E(5)*CSA(5) * D51(N5,:)); % pW5 w


% -B*v5 - A*w5 + v6 = 0
Gam(wC56(1),v6L) = 1; % v6

Gam(wC56(1),v5R) = -B; % v5
Gam(wC56(1),w5R) = -A; % w5

% B*pW5 - A*pV5 + pU6 = 0
Gam(wC56(2),u6L:u6R) = -E(6)*I(6)*D63(1,:); % pU6

Gam(wC56(2),u5R) = -B*E(5)*CSA(5) * kap5; % pW5 u
Gam(wC56(2),w5L:w5R) = B*E(5)*CSA(5) * D51(N5,:); % pW5 w

Gam(wC56(2),u5L:u5R) = Gam(wC56(2),u5L:u5R) + (3*A*E(5)*I(5) * kap5*mu5 * D52(N5,:)); % pV5 u
Gam(wC56(2),u5R) = Gam(wC56(2),u5R) - (A*E(5)*I(5) * (kap5^3)*(mu5^3)) - (2*A*G(5)*I(5) * (kap5^3)*mu5); % pV5 u
Gam(wC56(2),v5L:v5R) = -(A*E(5)*I(5) * ((3*(kap5^2)*(mu5^2)*D51(N5,:)) - D53(N5,:))) - (2*A*G(5)*I(5) * (kap5^2)*D51(N5,:)); % pV5 v
Gam(wC56(2),w5L:w5R) = Gam(wC56(2),w5L:w5R) + (2*A*E(5)*I(5) * (kap5^2)*mu5 * D51(N5,:)); % pV5 w
Gam(wC56(2),thet5L:thet5R) = -(A*E(5)*I(5) * kap5*D51(N5,:)) - (2*A*G(5)*I(5) * kap5*D51(N5,:)); % pV5 thetW

% B*thetW5 - A*thetV5 + thetU6 = 0
Gam(thetC56(1),v6L:v6R) = -D61(1,:); % thetU6

Gam(thetC56(1),thet5R) = B; % thetW5

Gam(thetC56(1),u5L:u5R) = -A*D51(N5,:); % thetV5 u
Gam(thetC56(1),v5R) = A*kap5*mu5; % thetV5 v
Gam(thetC56(1),w5R) = -A*kap5; % thetV5 w

% B*mW5 - A*mV5 + mU6 = 0
Gam(thetC56(2),v6L:v6R) = -E(6)*I(6)*D62(1,:); % mU6

Gam(thetC56(2),u5R) = 2*B*G(5)*I(5) * (kap5^2)*mu5; % mW5 u 
Gam(thetC56(2),v5L:v5R) = 2*B*G(5)*I(5) * kap5 * D51(N5,:); % mW5 v
Gam(thetC56(2),thet5L:thet5R) = 2*B*G(5)*I(5) * D51(N5,:); % mW5 thetW

Gam(thetC56(2),u5L:u5R) = Gam(thetC56(2),u5L:u5R) - (A*E(5)*I(5) * D52(N5,:)); % mV5 u
Gam(thetC56(2),u5R) = Gam(thetC56(2),u5R) + (A*E(5)*I(5) * (kap5^2)*(mu5^2)); % mV5 u
Gam(thetC56(2),v5L:v5R) = Gam(thetC56(2),v5L:v5R) + (2*A*E(5)*I(5)*kap5*mu5 * D51(N5,:)); % mV5 v
Gam(thetC56(2),w5L:w5R) = -A*E(5)*I(5)*kap5 * D51(N5,:); % mV5 w


%% CONTINUITY: SYSTEMS 6 -> 7 %%
% zero out rows to be replaced
Gam(uC67,:) = 0;
Gam(vC67,:) = 0;
Gam(wC67,:) = 0;
Gam(thetC67,:) = 0;

% zero for damping matrix too
R(uC67,:) = 0;
R(vC67,:) = 0;
R(wC67,:) = 0;
R(thetC67,:) = 0;

% u6 - u7 = 0
Gam(uC67(1),u6R) = 1;
Gam(uC67(1),u7L) = -1;
% thetV6 - thetW7 = 0
Gam(uC67(2),u6L:u6R) = D61(N6,:);
Gam(uC67(2),thet7L) = -1; % thetW7
% mV6 - mW7 = 0
Gam(uC67(3),u6L:u6R) = E(6)*I(6) * D62(N6,:); % mV6
Gam(uC67(3),thet7L:thet7R) = -2*G(7)*I(7) * D71(1,:); % mW7 
% pU6 - pU7 = 0
Gam(uC67(4),u6L:u6R) = -E(6)*I(6) * D63(N6,:); % pU6
Gam(uC67(4),u7L:u7R) = E(7)*I(7) * D73(1,:); % pU7

% v6 - w7 = 0
Gam(vC67(1),v6R) = 1;
Gam(vC67(1),w7L) = -1; 
% thetU6 - thetU7 = 0
Gam(vC67(2),v6L:v6R) = -D61(N6,:); % thetU6
Gam(vC67(2),v7L:v7R) = D71(1,:); % thetU7
% mU6 - mU7 = 0
Gam(vC67(3),v6L:v6R) = -E(6)*I(6) * D62(N6,:); % mU6
Gam(vC67(3),v7L:v7R) = E(7)*I(7) * D72(1,:); % mU7
% pV6 - pW7 = 0
Gam(vC67(4),v6L:v6R) = -E(6)*I(6) * D63(N6,:); % pV6
Gam(vC67(4),w7L:w7R) = -E(7)*CSA(7) * D71(1,:); % pW7

% w6 + v7 = 0
Gam(wC67(1),w6R) = 1;
Gam(wC67(1),v7L) = 1;
% pW6 + pV7 = 0
Gam(wC67(2),w6L:w6R) = E(6)*CSA(6) * D61(N6,:); % pW6
Gam(wC67(2),v7L:v7R) = -E(7)*I(7) * D73(1,:); % pV7

% thetW6 + thetV7 = 0
Gam(thetC67(1),thet6R) = 1; % thetW6
Gam(thetC67(1),u7L:u7R) = D71(1,:);
% mW6 + mV7 = 0
Gam(thetC67(2),thet6L:thet6R) = 2*G(6)*I(6) * D61(N6,:); % mW6
Gam(thetC67(2),u7L:u7R) = E(7)*I(7) * D72(1,:); % mV7


%% CONTINUITY: SYSTEMS 7 -> 8 %%
% zero out rows to be replaced
Gam(uC78,u7L:u8R) = 0;
Gam(vC78,v7L:v8R) = 0;
Gam(wC78,w7L:w8R) = 0;
Gam(thetC78,thet7L:thet8R) = 0;

% zero for damping matrix too
R(uC78,u7L:u8R) = 0;
R(vC78,v7L:v8R) = 0;
R(wC78,w7L:w8R) = 0;
R(thetC78,thet7L:thet8R) = 0;

% u
Gam(uC78(1),u7R) = 1;
Gam(uC78(1),u8L) = -1;

Gam(uC78(2),u7L:u7R) = D71(N7,:);
Gam(uC78(2),u8L:u8R) = -D81(1,:);

Gam(uC78(3),u7L:u7R) = E(7)*I(7) * D72(N7,:);
Gam(uC78(3),u8L:u8R) = E(8)*I(8) * -D82(1,:);

Gam(uC78(4),u7L:u7R) = E(7)*I(7) * D73(N7,:);
Gam(uC78(4),u8L:u8R) = E(8)*I(8) * -D83(1,:);

% v
Gam(vC78(1),v7R) = 1;
Gam(vC78(1),v8L) = -1;

Gam(vC78(2),v7L:v7R) = D71(N7,:);
Gam(vC78(2),v8L:v8R) = -D81(1,:);

Gam(vC78(3),v7L:v7R) = E(7)*I(7) * D72(N7,:);
Gam(vC78(3),v8L:v8R) = E(8)*I(8) * -D82(1,:);

Gam(vC78(4),v7L:v7R) = E(7)*I(7) * D73(N7,:);
Gam(vC78(4),v8L:v8R) = E(8)*I(8) * -D83(1,:);

% w
Gam(wC78(1),w7R) = 1;
Gam(wC78(1),w8L) = -1;

Gam(wC78(2),w7L:w7R) = E(7)*CSA(7) * D71(N7,:);
Gam(wC78(2),w8L:w8R) = E(8)*CSA(8) * -D81(1,:);

% theta_w
Gam(thetC78(1),thet7R) = 1;
Gam(thetC78(1),thet8L) = -1;

Gam(thetC78(2),thet7L:thet7R) = 2*G(7)*I(7) * D71(N7,:);
Gam(thetC78(2),thet8L:thet8R) = 2*G(8)*I(8) * -D81(1,:);


%% CONTINUITY: SYSTEMS 8 -> 9 %%
% zero out rows to be replaced
Gam(uC89,u8L:u9R) = 0;
Gam(vC89,v8L:v9R) = 0;
Gam(wC89,w8L:w9R) = 0;
Gam(thetC89,thet8L:thet9R) = 0;

% zero for damping matrix too
R(uC89,u8L:u9R) = 0;
R(vC89,v8L:v9R) = 0;
R(wC89,w8L:w9R) = 0;
R(thetC89,thet8L:thet9R) = 0;

% u
Gam(uC89(1),u8R) = 1;
Gam(uC89(1),u9L) = -1;

Gam(uC89(2),u8L:u8R) = D81(N8,:);
Gam(uC89(2),u9L:u9R) = -D91(1,:);

Gam(uC89(3),u8L:u8R) = E(8)*I(8) * D82(N8,:);
Gam(uC89(3),u9L:u9R) = E(9)*I(9) * -D92(1,:);

Gam(uC89(4),u8L:u8R) = E(8)*I(8) * D83(N8,:);
Gam(uC89(4),u9L:u9R) = E(9)*I(9) * -D93(1,:);

% v
Gam(vC89(1),v8R) = 1;
Gam(vC89(1),v9L) = -1;

Gam(vC89(2),v8L:v8R) = D81(N8,:);
Gam(vC89(2),v9L:v9R) = -D91(1,:);

Gam(vC89(3),v8L:v8R) = E(8)*I(8) * D82(N8,:);
Gam(vC89(3),v9L:v9R) = E(9)*I(9) * -D92(1,:);

Gam(vC89(4),v8L:v8R) = E(8)*I(8) * D83(N8,:);
Gam(vC89(4),v9L:v9R) = E(9)*I(9) * -D93(1,:);

% w
Gam(wC89(1),w8R) = 1;
Gam(wC89(1),w9L) = -1;

Gam(wC89(2),w8L:w8R) = E(8)*CSA(8) * D81(N8,:);
Gam(wC89(2),w9L:w9R) = E(9)*CSA(9) * -D91(1,:);

% theta_w
Gam(thetC89(1),thet8R) = 1;
Gam(thetC89(1),thet9L) = -1;

Gam(thetC89(2),thet8L:thet8R) = 2*G(8)*I(8) * D81(N8,:);
Gam(thetC89(2),thet9L:thet9R) = 2*G(9)*I(9) * -D91(1,:);


%% GROUP BOUNDARY AND CONTINUITY CONDITIONS %%
GamRowsMvd = [Gam(uB,:); Gam(vB,:); Gam(wB,:); Gam(thetB,:); ...
              Gam(uC,:); Gam(vC,:); Gam(wC,:); Gam(thetC,:); ...
              Gam(ui,:); Gam(vi,:); Gam(wi,:); Gam(theti,:)];

GamColsMvd = [GamRowsMvd(:,uB), GamRowsMvd(:,vB), GamRowsMvd(:,wB), GamRowsMvd(:,thetB), ...
              GamRowsMvd(:,uC), GamRowsMvd(:,vC), GamRowsMvd(:,wC), GamRowsMvd(:,thetC), ...
              GamRowsMvd(:,ui), GamRowsMvd(:,vi), GamRowsMvd(:,wi), GamRowsMvd(:,theti)];
          
RRowsMvd = [R(uB,:); R(vB,:); R(wB,:); R(thetB,:); ...
              R(uC,:); R(vC,:); R(wC,:); R(thetC,:); ...
              R(ui,:); R(vi,:); R(wi,:); R(theti,:)];

RColsMvd = [RRowsMvd(:,uB), RRowsMvd(:,vB), RRowsMvd(:,wB), RRowsMvd(:,thetB), ...
              RRowsMvd(:,uC), RRowsMvd(:,vC), RRowsMvd(:,wC), RRowsMvd(:,thetC), ...
              RRowsMvd(:,ui), RRowsMvd(:,vi), RRowsMvd(:,wi), RRowsMvd(:,theti)];


%% PARTITION TO SOLVE FOR INTERIOR %%
numBCs = 12; numCCs = 96;
numConds = numBCs + numCCs;

GamBB = GamColsMvd(1:numConds,1:numConds);
GamBI = GamColsMvd(1:numConds,(numConds+1):end);
GamIB = GamColsMvd((numConds+1):end,1:numConds);
GamII = GamColsMvd((numConds+1):end,(numConds+1):end);

% RBB, RBI are null matrices
RIB = RColsMvd((numConds+1):end,1:numConds);
RII = RColsMvd((numConds+1):end,(numConds+1):end);

Zw = (-GamIB * (GamBB \ GamBI)) + GamII;
Zr = (-RIB * (GamBB \ GamBI)) + RII;

Z = [zeros(size(Zw)), eye(size(Zw)); Zw, Zr];


%% EXTRACT EIGENVALUES & EIGENVECTORS %%
[P,Q] = eig(Z);