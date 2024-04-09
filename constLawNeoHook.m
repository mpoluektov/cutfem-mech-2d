function [ P, d_PT_d_F, W, d_AderB_d_F ] = constLawNeoHook( F, K, G, g, Am, Bm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

iG = diag( [ 1/g 1/g 1 ] );
detG = g^2;
Fe = F * iG;

Je = det(Fe);
FeT = Fe.';
Be = Fe * FeT;
Bei = (Je^(-2/3)) * Be;
Beid = Bei - (1/3) * tensTrace3(Bei) * eye(3);
tauh = K * ( Je - 1 ) * eye(3);
taud = G * Beid;
tau = tauh + taud;
iFe = tensInv3(Fe);
iFeT = iFe.';
Pe = tau * iFeT;

P = detG * Pe * iG.';

W = detG * ( K * ( Je - 1 - log(Je) ) + (G/2) * ( tensTrace3(Bei) - 3 ) );

%% derivatives

C3 = zeros(9);
C3(1:3,1:3) = eye(3);
C3(4:6,7:9) = eye(3);
C3(7:9,4:6) = eye(3);

iG4I = tens2mlt4I(iG);
Fe4I = tens2mlt4I(Fe);
Pe4I = tens2mlt4I(Pe);
iFe4I = tens2mlt4I(iFe);
vFe = tens2vec(Fe);
vFeT = tens2vec(FeT);
viFe = tens2vec(iFe);
viFeT = tens2vec(iFeT);
GJe23 = G * (Je^(-2/3));

d_PeT_d_Fe = K * Je * viFe * viFeT.' + GJe23 * ( C3 - (2/3) * vFeT * viFeT.' + iFe4I * Fe4I + ...
             (2/9) * tensTrace3( Fe * FeT ) * viFe * viFeT.' - (2/3) * viFe * vFe.' ) - iFe4I * Pe4I;

d_PT_d_F = detG * iG4I * C3 * d_PeT_d_Fe * iG4I.' * C3;

%% second derivatives

%. calculate d(Am:d(PT)/d(F):Bm)/d(F)

if ( nargin == 6 )

    A = detG * Am * iG;
    B = iG.' * Bm;

    iFeTAT = iFeT * A.';
    BiFeT = B * iFeT;
    trAiFe = tensTrace3( iFeTAT );
    trAFeT = tensTrace3( A * FeT );
    triFeTB = tensTrace3( BiFeT );
    trFeB = tensTrace3( Fe * B );
    trAB = tensTrace3( A * B );
    trBe = tensTrace3( Be );
    AiFeBT = A * BiFeT.';
    AiFeBTFeT = AiFeBT * FeT;
    trAiFeBTFeT = tensTrace3( AiFeBTFeT );
    vAiFeBT = tens2vec( AiFeBT );
    vAiFeBTder = d_PeT_d_Fe.' * C3 * vAiFeBT;
    AiFeBTder = vec2tens( vAiFeBTder );
    iFeTATiFeT = iFeTAT * iFeT;
    iFeTBiFeT = iFeT * BiFeT;

    d_AderB_d_Fe = K * Je * trAiFe * triFeTB * iFeT - K * Je * triFeTB * iFeTATiFeT - K * Je * trAiFe * iFeTBiFeT + ...
                   (-2/3) * GJe23 * trAB * iFeT + ...
                   (4/9) * GJe23 * trAFeT * triFeTB * iFeT + (-2/3) * GJe23 * triFeTB * A + (2/3) * GJe23 * trAFeT * iFeTBiFeT + ...
                   (-2/3) * GJe23 * trAiFeBTFeT * iFeT - GJe23 * iFeTAT * Fe * BiFeT + GJe23 * AiFeBT + ...
                   (-4/27) * GJe23 * trBe * trAiFe * triFeTB * iFeT + (4/9) * GJe23 * trAiFe * triFeTB * Fe + ...
                   (-2/9) * GJe23 * trBe * triFeTB * iFeTATiFeT + (-2/9) * GJe23 * trBe * trAiFe * iFeTBiFeT + ...
                   (4/9) * GJe23 * trAiFe * trFeB * iFeT + (2/3) * GJe23 * trFeB * iFeTATiFeT + (-2/3) * GJe23 * trAiFe * B.' + ...
                   iFeTAT * Pe * BiFeT - AiFeBTder;

    d_AderB_d_F = tens2vec(d_AderB_d_Fe * iG.').';
    
end

end

