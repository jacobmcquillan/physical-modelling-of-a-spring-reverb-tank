%% LOSSY SPRING REVERB MODEL %%
% JACOB MCQUILLAN
% SARC, QUEEN'S UNIVERSITY BELFAST

% REQUIREMENTS: install chebfun: chebfun.org
% codes uses diffmat.m to generate differentiation matrices and chebpts.m
% to generate Chebyshev grids

% set model parameters in reverbTankParams.m
% eigenvalues extracted in tankModelEig.m
% modal algorithm runs in this code

clear;


%% SET-UP %%
% choose boundary conditions
f_BCs = "clamped"; % clamped or pinned

f_velOut = 0; % 1: velocity as output (const I source)
              % 0: displacement out (const V source)
              % just re-run from mode amplitude computation for quick comp

% maximum frequency of modes to retain
flim = 12e3;

% set grid points
N1 = 30; % all wires set the same for convenience
N2 = 800; % helical spring (800 reasonable for Olson)
N = [N1 N1 N1 N1 N2 N1 N1 N1 N1]; % consolidate

% damping parameters
sigA = 0.4;
sigB = 4;
sigC = 5;
eta = 180;


%% DISCRETISE MODEL & PERFORM EIGENVALUE ANALYSIS %%
[P,Q,s,Zw,Zr] = tankModelEig(f_BCs,N,eta,sigA,sigB,sigC);
% P: columns of eigenvectors corresponding to ...
% Q: diagonal eigenvalue matrix (eigenvalues appear as complex conjugate
% pairs)
% s: grid points of the spatial domains concatenated
% Zw, Zr: state-space system matrices, Z = [0, I; Zw, Zr]

% load material and geometry (convenient to have in workspace)
load SpringReverbParameters


%% COMPLEX DIAGONAL FORM TO REAL BLOCK DIAGONAL FORM %%
[Ptild,Qtild] = cdf2rdf(P,Q);
% lam = Re(lam) +/- imag(lam) -> [Re(lam), Imag(lam); -Imag(lam), Re(lam)]
% Ptild satisfies Z*Ptild = Ptild*Qtild


%% SPATIAL GRIDS %%
s1 = s(1:N1);
s2 = s(N1+1:2*N1);
s3 = s((2*N1)+1:3*N1);
s4 = s((3*N1)+1:(3*N1)+N1);
s5 = s(((3*N1)+N1+1):(3*N1)+N1+N2);
s6 = s(((3*N1)+N1+N2+1):(3*N1)+(2*N1)+N2);
s7 = s((3*N1)+(2*N1)+N2+1:(4*N1)+(2*N1)+N2);
s8 = s((4*N1)+(2*N1)+N2+1:(5*N1)+(2*N1)+N2);
s9 = s((5*N1)+(2*N1)+N2+1:(6*N1)+(2*N1)+N2);

% interior grids for each polarisation
sUi = [s1(3:end-2); L(1)+s2(3:end-2); sum(L(1:2))+s3(3:end-2); ...
       sum(L(1:3))+s4(3:end-2); sum(L(1:4))+s5(3:end-2); ...
       sum(L(1:5))+s6(3:end-2); sum(L(1:6))+s7(3:end-2); ...
       sum(L(1:7))+s8(3:end-2); sum(L(1:8))+s9(3:end-2);];
sVi = sUi;
sWi = [s1(2:end-1); L(1)+s2(2:end-1); sum(L(1:2))+s3(2:end-1); ...
       sum(L(1:3))+s4(2:end-1); sum(L(1:4))+s5(2:end-1); ...
       sum(L(1:5))+s6(2:end-1); sum(L(1:6))+s7(2:end-1); ...
       sum(L(1:7))+s8(2:end-1); sum(L(1:8))+s9(2:end-1)];
sTheti = sWi;

% concatenated interior grids
si = [sUi; sum(L)+sVi; (2*sum(L))+sWi; (3*sum(L))+sTheti];


%% DEFINE EXCITATION AND PICK-UP %%
% excite theta_w at left bead and pick-up theta_w at right bead
psiE = zeros(length(si),1); psiP = psiE;

% need to distribute evenly across the Chebyshev grid
wght = sin((pi*(0:N1-1))/(N1-1));
wghti = wght(2:(end-1));

% excitation distribution
psiE(si>((3*sum(L))+L(1))) = 1 / L(2);
psiE(si>((3*sum(L))+sum(L(1:2)))) = 0;

psiE(psiE~=0) = psiE(psiE~=0) .* wghti'; % weighting

psiE = [zeros(size(psiE)); psiE]; % state-space form

% pick-up distribution
psiP(si>((3*sum(L))+sum(L(1:7)))) = 1 / L(8);
psiP(si>((3*sum(L))+sum(L(1:8)))) = 0;

psiP(psiP~=0) = psiP(psiP~=0) .* wghti'; % weighting


%% MODE AMPLITUDES %%
% input modal amplitudes
cI = Ptild \ psiE;

% output modal amplitudes
if f_velOut
    cO = [zeros(1,length(psiP)), psiP'] * Ptild; % velocity output
else
    cO = [psiP', zeros(1,length(psiP))] * Ptild; % displacement output
end


%% SELECT FREQUENCY RANGE SET ON LINE 19 %%
dec = -diag(Qtild); % decay rates
om = zeros(length(Qtild),1); % frequencies [rad]

% this looks a bit convoluted but it makes the next section easier
for n=1:length(om)
    if mod(n,2) ~= 0 % n = 1,3,5...
        om(n)  = Qtild(n,n+1);
    else % n = 2,4,6...
        om(n)  = Qtild(n,n-1);
    end
end

f = om ./ (2*pi); % frequencies [Hz]
mat = [om,dec,cI,cO']; % consolidate

% cut anything outside of frequency range defined at start
mataud = mat(abs(om)<=2*pi*flim,:);


%% MODE PARAMETER EXTRACTION %%
% modes appear in 2x2 blocks, but overdamped modes appear in a 1x1 block so
% we need to pad those into 2x2 blocks to makes things easier

m = 0;
for n=1:length(mataud)
    n = n+m; % an extra increment to skip rows we just added
    if mataud(n,1) == 0 % overdamped mode (only this simple bc of how  prev. section is set up)
        mataud = [mataud(1:n,:); mataud(n,:); mataud(n+1:end,:)];
        m = m+1; % skip the row we just added next time
    end
end

omaud = mataud(1:2:end,1); % frequencies
decaud = mataud(1:2:end,2); % decay rates
cIp = mataud(1:2:end,3); cIq = mataud(2:2:end,3); % input amplitudes
cOp = mataud(1:2:end,4); cOq = mataud(2:2:end,4); % output amplitudes

faud = omaud/(2*pi);


%% TIME-DOMAIN SETUP %%
Fs = 48e3; % sampling frequency [Hz]
dt = 1/Fs; % corresponding time step [s]
dur = 3; % signal duration [s]
Ns = dur*Fs; % no. of samples
t = (0:(Ns-1))*dt; % time vector


%% PARAMETERS FOR IMPULSE INVARIANT METHOD %%
% update equation coefficients
b0p = dt*cIp;
b1p = dt*exp(-decaud*dt) .* (cIq.*sin(omaud*dt) - cIp.*cos(omaud*dt));

b0q = dt*cIq;
b1q = dt*exp(-decaud*dt) .* (-cIp.*sin(omaud*dt) - cIq.*cos(omaud*dt));

a1 = -2*exp(-decaud*dt).*cos(omaud*dt);
a2 = exp(-2*decaud*dt);

% only keep dissipative modes
b0p(decaud<0) = 0; b1p(decaud<0) = 0;
b0q(decaud<0) = 0; b1q(decaud<0) = 0;
a1(decaud<0) = 0; a2(decaud<0) = 0;


%% MODAL ALGORITHM %%
Fe = [1; zeros(Ns-1,1)]; % input

p = zeros(length(b0p),1); pp = p; pm = p;
q = p; qp = p; qm = p;
xp = 0; x = 0;

Ns_loop = Ns; % sometimes convenient to only compute some
    
Fp = zeros(Ns_loop,1); % output

for n=1:Ns_loop
    xp = Fe(n);

    % update equations
    pp = b0p.*xp + b1p.*x - a1.*p - a2.*pm;
    qp = b0q.*xp + b1q.*x - a1.*q - a2.*qm;
    
    % compute output
    Fp(n) = (cOp' * pp) + (cOq' * qp);
    
    % update for next time step
    pm = p; p = pp;
    qm = q; q = qp;
    x = xp;
end


%% PLOTTING %%
figure(1); clf;
plot(t,Fp);