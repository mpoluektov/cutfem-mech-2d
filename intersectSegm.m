function [ S ] = intersectSegm(c1,c2,p1,p2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A = [ c2(1)-c1(1)  p1(1)-p2(1) ;
      c2(2)-c1(2)  p1(2)-p2(2) ];
B = [ p1(1)-c1(1) ;
      p1(2)-c1(2) ];
S = [ A(2,2)*B(1) - A(1,2)*B(2) ;
      A(1,1)*B(2) - A(2,1)*B(1) ] / ( A(1,1)*A(2,2) - A(1,2)*A(2,1) );

end

