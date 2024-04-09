function [ P, d_PT_d_F, W, Cp, d_AderB_d_F ] = constLawVisc( F, Cp_prev, K, G, eta, g, dt, Am, Bm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% chemical part

iG = diag( [ 1/g 1/g 1 ] );
detG = g^2;
Fr = F * iG;

%% plastic flow

Je = det(Fr);
Cr = Fr.' * Fr;
Cp_inc = Cr * (G*dt/eta) * Je^(-2/3) + Cp_prev;
Cp = Cp_inc * det(Cp_inc)^(-1/3);
iCp = tensInv3(Cp);
Be = Fr * iCp * Fr.';

%% elasticity

Bei = Be * Je^(-2/3);
Beid = Bei - (1/3) * tensTrace3(Bei) * eye(3);
tauh = K * ( Je - 1 ) * eye(3);
taud = G * Beid;
tau = tauh + taud;

%% Piola-Kirchhoff and energy

iFr = tensInv3(Fr);
Pr = tau * iFr.';
Wr = K * ( Je - 1 - log(Je) ) + (G/2) * ( tensTrace3(Bei) - 3 );
P = detG * Pr * iG.';
W = detG * Wr;

%% derivatives

%. tensor product of 2nd-order identity tensors
C1 = zeros(9);
C1(1:3,1:3) = ones(3);

%. right-transpose of 4th-order identity tensor
C2 = eye(9);

%. 4th-order identity tensor
C3 = zeros(9);
C3(1:3,1:3) = eye(3);
C3(4:6,7:9) = eye(3);
C3(7:9,4:6) = eye(3);

%. d(tau)/d(Je)
taus = ( tauh + K * eye(3) ) / Je + taud * (-2/3) / Je;

iG4I = tens2mlt4I(iG);
Fr4I = tens2mlt4I(Fr);
iFr4I = tens2mlt4I(iFr);
iCp4I = tens2mlt4I(iCp);
tau4I = tens2mlt4I(tau);
vCr = tens2vec(Cr);
viFr = tens2vec(iFr);
vCp = tens2vec(Cp);
viCp = tens2vec(iCp);
vtaus = tens2vec(taus);

d_Cr_d_Fr = ( C2 + C3 ) * Fr4I.';
d_Cp_inc_d_Fr = ( d_Cr_d_Fr + (-2/3) * vCr * viFr.' * C3 ) * (G*dt/eta) * Je^(-2/3);
d_Cp_d_Fr = ( C2 + (-1/3) * vCp * viCp.' ) * d_Cp_inc_d_Fr * det(Cp_inc)^(-1/3);
d_iCp_d_Cp = -iCp4I * iCp4I;
d_iCp_d_Fr = d_iCp_d_Cp * d_Cp_d_Fr;
d_Be_d_Fr = ( C2 + C3 ) * Fr4I * C3 * iCp4I + Fr4I * Fr4I * d_iCp_d_Fr;
d_tau_d_Be = G * Je^(-2/3) * ( C2 - (1/3) * C1 );
d_tau_d_Fr = d_tau_d_Be * d_Be_d_Fr + vtaus * Je * viFr.' * C3;
d_iFr_d_Fr = -iFr4I * C3 * iFr4I.' * C3;
d_Pr_d_Fr = tau4I * d_iFr_d_Fr + C3 * iFr4I * d_tau_d_Fr;
d_PT_d_F = detG * iG4I * d_Pr_d_Fr * iG4I.' * C3;

%% second derivatives

%. calculate d(Am:d(PT)/d(F):Bm)/d(F)

if ( nargin == 9 )

    vAmT = tens2vec(Am.');
    vBmT = tens2vec(Bm.');
    R = vAmT.' * d_PT_d_F * vBmT ;

    SDLT = 1e-8;
    d_AderB_d_F_r = zeros(3,3);
    for ii=1:9
        F_c = F;
        F_c(ii) = F_c(ii) + SDLT;
        [ dum, d_PT_d_F_c ] = constLaw_visc( F_c, Cp_prev, K, G, eta, g, dt );
        R_c = vAmT.' * d_PT_d_F_c * vBmT;
        d_AderB_d_F_r(ii) = ( R_c - R ) / SDLT;
    end
    d_AderB_d_F = tens2vec(d_AderB_d_F_r).';

end

end

