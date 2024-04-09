function [ S ] = intersectCurve( origPoints, origConn, p1, p2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nsegm = size( origConn, 1 );

S = zeros(2,1);

dist_in = Inf;
dist_out = Inf;
foundInside = 0;
for jj=1:Nsegm
    c1 = origPoints(origConn(jj,1),:);
    c2 = origPoints(origConn(jj,2),:);
    st = intersectSegm(c1,c2,p1,p2);
    dist_test = ( st(1) - 0.5 )^2 + ( st(2) - 0.5 )^2;
    doReassign = 0;
    if ( st(1) >= 0 )&&( st(1) <= 1 )&&( st(2) >= 0 )&&( st(2) <= 1 )
        if ( dist_test < dist_in )
            dist_in = dist_test;
            foundInside = 1;
            doReassign = 1;
        end
    elseif ( dist_test < dist_out )&&( foundInside == 0 )
        dist_out = dist_test;
        doReassign = 1;
    end
    if ( doReassign == 1 )
        S = st;
    end
end

end

