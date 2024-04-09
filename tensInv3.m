function [ iA ] = tensInv3( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dA = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1);
iA = [ A(2,2)*A(3,3)-A(2,3)*A(3,2)  A(1,3)*A(3,2)-A(1,2)*A(3,3)  A(1,2)*A(2,3)-A(1,3)*A(2,2) ;
       A(2,3)*A(3,1)-A(2,1)*A(3,3)  A(1,1)*A(3,3)-A(1,3)*A(3,1)  A(1,3)*A(2,1)-A(1,1)*A(2,3) ;
       A(2,1)*A(3,2)-A(2,2)*A(3,1)  A(1,2)*A(3,1)-A(1,1)*A(3,2)  A(1,1)*A(2,2)-A(1,2)*A(2,1) ]/dA;

end

