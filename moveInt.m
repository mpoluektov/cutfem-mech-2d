function [ intPoints_new ] = moveInt( conc, ceq, xi, matPar, intParam, dt )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

coefVel = matPar.coefVel;

fracFaces = intParam.fracFaces;
intNodes = intParam.intNodes;
intPoints = intParam.intPoints;
intNorms = intParam.intNorms;

NintP = size(intPoints,1);

%. velocities
if ( matPar.isChem == 1 )
    %. chemo-mechanics
    conc_int = zeros(NintP,1);
    for ii=1:NintP
        n1 = intNodes(ii,1);
        n2 = intNodes(ii,2);
        w1 = fracFaces(ii);
        w2 = 1 - fracFaces(ii);
        cn1 = conc(n1);
        cn2 = conc(n2);
        conc_int(ii) = cn1*w1 + cn2*w2;
    end
    vel_int = coefVel * (conc_int - ceq);
else
    %. phase transition
    vel_int = coefVel * xi;
end

%. new positions of the points
intPoints_new = intPoints + dt * intNorms .* repmat( vel_int, 1, 2 );

end

