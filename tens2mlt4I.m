function [ A4I ] = tens2mlt4I( A )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%. evaluate: A \cdot {^4}I

vA = tens2vec(A);

A4I = [ vA(1)  0  0  vA(4)  vA(5)  0  0  0  0 ;
        0  vA(2)  0  0  0  vA(6)  vA(7)  0  0 ;
        0  0  vA(3)  0  0  0  0  vA(8)  vA(9) ;
        0  vA(4)  0  0  0  vA(5)  vA(1)  0  0 ;
        0  0  vA(5)  0  0  0  0  vA(1)  vA(4) ;
        0  0  vA(6)  0  0  0  0  vA(7)  vA(2) ;
        vA(7)  0  0  vA(2)  vA(6)  0  0  0  0 ;
        vA(8)  0  0  vA(9)  vA(3)  0  0  0  0 ;
        0  vA(9)  0  0  0  vA(3)  vA(8)  0  0 ];

end

