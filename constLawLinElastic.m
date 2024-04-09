function [ P, d_PT_d_F, W ] = constLawLinElastic( F, K, G, g )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

I = eye(3);
epsch = diag( [ g g 1 ] ) - I;

eps = ( F + F.' - 2*I )/2 - epsch;
sig = (K-(2/3)*G)*tensTrace3(eps)*I + 2*G*eps;

P = sig;
W = tensTrace3(sig*eps)/2;

%% derivatives

C1 = zeros(9);
C1(1:3,1:3) = ones(3);
C2 = eye(9);
C3 = zeros(9);
C3(1:3,1:3) = eye(3);
C3(4:6,7:9) = eye(3);
C3(7:9,4:6) = eye(3);

d_sig_d_Fe = (K-(2/3)*G)*C1 + G*(C2+C3);
d_PT_d_F = d_sig_d_Fe;

end

