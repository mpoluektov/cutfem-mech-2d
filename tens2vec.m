function [ B ] = tens2vec( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

B = [ A(1,1) A(2,2) A(3,3) A(1,2) A(1,3) A(2,3) A(2,1) A(3,1) A(3,2) ].';

end

