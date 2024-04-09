function [ P, d_PT_d_F, W, Cp, d_AderB_d_F ] = constLaw( F, Cp_prev, K, G, eta, g, dt, Am, Bm )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ( nargin == 9 )
%     [ P, d_PT_d_F, W, Cp, d_AderB_d_F ] = constLawVisc( F, Cp_prev, K, G, eta, g, dt, Am, Bm );
%     [ P, d_PT_d_F, W ] = constLawLinElastic( F, K, G, g, Am, Bm );
%     d_AderB_d_F = zeros(1,9);
    [ P, d_PT_d_F, W, d_AderB_d_F ] = constLawNeoHook( F, K, G, g, Am, Bm );
    Cp = Cp_prev;
else
%     [ P, d_PT_d_F, W, Cp ] = constLawVisc( F, Cp_prev, K, G, eta, g, dt );
%     [ P, d_PT_d_F, W ] = constLawLinElastic( F, K, G, g );
    [ P, d_PT_d_F, W ] = constLawNeoHook( F, K, G, g );
    Cp = Cp_prev;
end

end

