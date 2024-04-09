function [ intParam ] = calcIntLevelSet( origPoints, origConn, geomPar, PTOL )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dh = geomPar.dh;

windThresh = geomPar.windThresh;

%% grid

gXs = ( 0:1:(Nx-1) )*dh;
gYs = ( 0:1:(Ny-1) )*dh;
gX = repmat(gXs,1,Ny);
gYa = repmat(gYs,1,Nx);
gYb = reshape(gYa,Ny,[]);
gYt = gYb.';
gY = gYt(:).';
coordGrid = [ gX; gY; ];

%% calculate winding numbers

Nsegm = size( origConn, 1 );

gridWind = zeros( Nx*Ny, 1 );
for ii=1:Nx*Ny
    for jj=1:Nsegm
        c1 = origPoints(origConn(jj,1),:);
        c2 = origPoints(origConn(jj,2),:);
        p1 = coordGrid(:,ii).';
        v1 = c1-p1;
        v2 = c2-p1;
        crProd = v1(1)*v2(2) - v1(2)*v2(1);
        dtProd = v1(1)*v2(1) + v1(2)*v2(2);
        thet = atan2(crProd,dtProd);
        gridWind(ii) = gridWind(ii) + thet;
    end
end
gridWindSc = gridWind/(2*pi);

%% assign element types

elTypes = zeros( 2*(Nx-1)*(Ny-1), 1 );
for ii=1:((Nx-1)*(Ny-1))

    ind_bl = ii + floor((ii-0.5)/(Nx-1));
    ind_br = ii + 1 + floor((ii-0.5)/(Nx-1));
    ind_tl = ii + Nx + floor((ii-0.5)/(Nx-1));
    ind_tr = ii + Nx + 1 + floor((ii-0.5)/(Nx-1));
    
    elem = [ ind_tl  ind_bl  ind_br ;
             ind_br  ind_tr  ind_tl ];
    
    for jj=1:2
        if ( gridWindSc(elem(jj,1)) <= windThresh )&&( gridWindSc(elem(jj,2)) <= windThresh )&&( gridWindSc(elem(jj,3)) <= windThresh )
            %. transf
            elTypes( 2*ii+jj-2 ) = 1;
        elseif ( gridWindSc(elem(jj,1)) >= windThresh )&&( gridWindSc(elem(jj,2)) >= windThresh )&&( gridWindSc(elem(jj,3)) >= windThresh )
            %. untransf
            elTypes( 2*ii+jj-2 ) = 2;
        else
            %. cut
            elTypes( 2*ii+jj-2 ) = 3;
        end
    end
end

%% find cut edges

intElem = find( elTypes == 3 );

NintE = size(intElem,1);
cutTypes = zeros(NintE,2);
cutInds = zeros(NintE,2);
for ii=1:NintE

    %. elem. idx
    eind = intElem(ii);
    
    %. square idx
    sqind = ceil((eind-0.5)/2);
    
    %. nodes that belong to the square
    ind_bl = sqind + floor((sqind-0.5)/(Nx-1));
    ind_br = sqind + 1 + floor((sqind-0.5)/(Nx-1));
    ind_tl = sqind + Nx + floor((sqind-0.5)/(Nx-1));
    ind_tr = sqind + Nx + 1 + floor((sqind-0.5)/(Nx-1));

    %. edges
    ind_L = ind_bl;
    ind_R = ind_br;
    ind_B = sqind;
    ind_T = sqind+Nx-1;
    ind_D = sqind;
    
    if ( mod(eind,2)==0 )
        %. top triangle
        if ( gridWindSc(ind_tl) <= windThresh )&&( gridWindSc(ind_tr) > windThresh )&&( gridWindSc(ind_br) < windThresh )
            %. T R
            cutTypes(ii,:) = [2 1];
            cutInds(ii,:) = [ind_T ind_R];
        elseif ( gridWindSc(ind_tl) < windThresh )&&( gridWindSc(ind_tr) <= windThresh )&&( gridWindSc(ind_br) > windThresh )
            %. R D
            cutTypes(ii,:) = [1 3];
            cutInds(ii,:) = [ind_R ind_D];
        elseif ( gridWindSc(ind_tl) < windThresh )&&( gridWindSc(ind_tr) > windThresh )&&( gridWindSc(ind_br) >= windThresh )
            %. T D
            cutTypes(ii,:) = [2 3];
            cutInds(ii,:) = [ind_T ind_D];
        elseif ( gridWindSc(ind_tl) >= windThresh )&&( gridWindSc(ind_tr) < windThresh )&&( gridWindSc(ind_br) > windThresh )
            %. R T
            cutTypes(ii,:) = [1 2];
            cutInds(ii,:) = [ind_R ind_T];
        elseif ( gridWindSc(ind_tl) > windThresh )&&( gridWindSc(ind_tr) >= windThresh )&&( gridWindSc(ind_br) < windThresh )
            %. D R
            cutTypes(ii,:) = [3 1];
            cutInds(ii,:) = [ind_D ind_R];
        elseif ( gridWindSc(ind_tl) > windThresh )&&( gridWindSc(ind_tr) < windThresh )&&( gridWindSc(ind_br) <= windThresh )
            %. D T
            cutTypes(ii,:) = [3 2];
            cutInds(ii,:) = [ind_D ind_T];
        end
    else
        %. bottom triangle
        if ( gridWindSc(ind_tl) <= windThresh )&&( gridWindSc(ind_bl) > windThresh )&&( gridWindSc(ind_br) < windThresh )
            %. B L
            cutTypes(ii,:) = [2 1];
            cutInds(ii,:) = [ind_B ind_L];
        elseif ( gridWindSc(ind_tl) < windThresh )&&( gridWindSc(ind_bl) <= windThresh )&&( gridWindSc(ind_br) > windThresh )
            %. D B
            cutTypes(ii,:) = [3 2];
            cutInds(ii,:) = [ind_D ind_B];
        elseif ( gridWindSc(ind_tl) < windThresh )&&( gridWindSc(ind_bl) > windThresh )&&( gridWindSc(ind_br) >= windThresh )
            %. D L
            cutTypes(ii,:) = [3 1];
            cutInds(ii,:) = [ind_D ind_L];
        elseif ( gridWindSc(ind_tl) >= windThresh )&&( gridWindSc(ind_bl) < windThresh )&&( gridWindSc(ind_br) > windThresh )
            %. L B
            cutTypes(ii,:) = [1 2];
            cutInds(ii,:) = [ind_L ind_B];
        elseif ( gridWindSc(ind_tl) > windThresh )&&( gridWindSc(ind_bl) >= windThresh )&&( gridWindSc(ind_br) < windThresh )
            %. B D
            cutTypes(ii,:) = [2 3];
            cutInds(ii,:) = [ind_B ind_D];
        elseif ( gridWindSc(ind_tl) > windThresh )&&( gridWindSc(ind_bl) < windThresh )&&( gridWindSc(ind_br) <= windThresh )
            %. L D
            cutTypes(ii,:) = [1 3];
            cutInds(ii,:) = [ind_L ind_D];
        end
    end
end

%% find intersection points

cutTypes_r = cutTypes(:);
cutInds_r = cutInds(:);
cutFacesVert = unique( cutInds_r(cutTypes_r==1) );
cutFacesHorz = unique( cutInds_r(cutTypes_r==2) );
cutFacesDiag = unique( cutInds_r(cutTypes_r==3) );

NVert = size(cutFacesVert,1);
NHorz = size(cutFacesHorz,1);
NDiag = size(cutFacesDiag,1);

fracFacesVert = zeros(NVert,1);
fracFacesHorz = zeros(NHorz,1);
fracFacesDiag = zeros(NDiag,1);

intPointsVert = zeros(NVert,2);
intPointsHorz = zeros(NHorz,2);
intPointsDiag = zeros(NDiag,2);

intNodesVert = zeros(NVert,2);
intNodesHorz = zeros(NHorz,2);
intNodesDiag = zeros(NDiag,2);

for ii=1:NVert
    indF = cutFacesVert(ii);
    indN = [ indF indF+Nx ];
    p1 = coordGrid(:,indN(1)).';
    p2 = coordGrid(:,indN(2)).';
    st = intersectCurve( origPoints, origConn, p1, p2 );
    fracFacesVert(ii) = 1-st(2);
    intNodesVert(ii,:) = indN;
    intPointsVert(ii,:) = p1 + st(2)*( p2 - p1 );
end

for ii=1:NHorz
    indF = cutFacesHorz(ii);
    indN = [ indF indF+1 ] + floor((indF-0.5)/(Nx-1));
    p1 = coordGrid(:,indN(1)).';
    p2 = coordGrid(:,indN(2)).';
    st = intersectCurve( origPoints, origConn, p1, p2 );
    fracFacesHorz(ii) = 1-st(2);
    intNodesHorz(ii,:) = indN;
    intPointsHorz(ii,:) = p1 + st(2)*( p2 - p1 );
end

for ii=1:NDiag
    indF = cutFacesDiag(ii);
    indN = [ indF+Nx indF+1 ] + floor((indF-0.5)/(Nx-1));
    p1 = coordGrid(:,indN(1)).';
    p2 = coordGrid(:,indN(2)).';
    st = intersectCurve( origPoints, origConn, p1, p2 );
    fracFacesDiag(ii) = 1-st(2);
    intNodesDiag(ii,:) = indN;
    intPointsDiag(ii,:) = p1 + st(2)*( p2 - p1 );
end

%% link cut elements and intersection points

fracFaces = [ fracFacesVert; fracFacesHorz; fracFacesDiag; ];
intNodes = [ intNodesVert; intNodesHorz; intNodesDiag; ];
intPoints = [ intPointsVert; intPointsHorz; intPointsDiag; ];

intConn = zeros(NintE,2);
for ii=1:NintE
    for jj=1:2
        if ( cutTypes(ii,jj) == 1 )
            locInd = find( cutFacesVert == cutInds(ii,jj), 1, 'first' );
            intConn(ii,jj) = locInd;
        elseif ( cutTypes(ii,jj) == 2 )
            locInd = find( cutFacesHorz == cutInds(ii,jj), 1, 'first' );
            intConn(ii,jj) = locInd + NVert;
        elseif ( cutTypes(ii,jj) == 3 )
            locInd = find( cutFacesDiag == cutInds(ii,jj), 1, 'first' );
            intConn(ii,jj) = locInd + NVert + NHorz;
        end
    end
end

%% element normals and volume fractions

cutLen = zeros(NintE,1);
elemNorms = zeros(NintE,2);
fracElem = zeros(NintE,1);
for ii=1:NintE

    intPointsA = intPoints(intConn(ii,1),:);
    intPointsB = intPoints(intConn(ii,2),:);
    secVec = intPointsB - intPointsA;
    secLen = norm(secVec);
    if ( secLen > 0 )
        secVecNrm = secVec/secLen;
        secNorm = [ -secVecNrm(2) secVecNrm(1) ];
    else
        secNorm = [ 1 0 ];
    end

    cutLen(ii) = secLen;
    elemNorms(ii,:) = secNorm;

    intNodesA = intNodes(intConn(ii,1),:);
    intNodesB = intNodes(intConn(ii,2),:);
    [comNode, iA, iB] = intersect(intNodesA,intNodesB);
    if ( iA == 1 )
        fracA = 1 - fracFaces(intConn(ii,1));
    else
        fracA = fracFaces(intConn(ii,1));
    end
    if ( iB == 1 )
        fracB = 1 - fracFaces(intConn(ii,2));
    else
        fracB = fracFaces(intConn(ii,2));
    end
    if ( cutTypes(ii,1) == 3 )
        fracA = fracA * sqrt(2);
    end
    if ( cutTypes(ii,2) == 3 )
        fracB = fracB * sqrt(2);
    end
    fracS = secLen/dh;
    areaTriang = (1/4)*sqrt(abs((fracA+fracB+fracS)*(-fracA+fracB+fracS)*(fracA-fracB+fracS)*(fracA+fracB-fracS)));
    if ( gridWindSc(comNode) < windThresh )
        fracElem(ii) = 1 - 2*areaTriang;
    else
        fracElem(ii) = 2*areaTriang;
    end
end

%% point normals

NintP = size(intPoints,1);
intNorms = zeros(NintP,2);
for ii=1:NintP
    curPoint = intPoints(ii,:);
    pDiff = intPoints - repmat(curPoint,NintP,1);
    pDist = sqrt( pDiff(:,1).^2 + pDiff(:,2).^2 );
    indC = find( pDist < PTOL );
    indA = find( ismember( intConn(:,1), indC ) | ismember( intConn(:,2), indC ) );
    currNorm = sum( elemNorms(indA,:) .* repmat(cutLen(indA,:),1,2), 1 );
    intNorms(ii,:) = currNorm/norm(currNorm);
end

%% excluded elements

intNodesList = unique(intNodes(:));
exclNodes_U = find( gridWindSc < windThresh );
indrem = ismember( exclNodes_U, intNodesList );
intNodesUnk_U = exclNodes_U(indrem);
exclNodes_U(indrem) = [];
exclNodes_T = find( gridWindSc > windThresh );
indrem = ismember( exclNodes_T, intNodesList );
intNodesUnk_T = exclNodes_T(indrem);
exclNodes_T(indrem) = [];
intNodesUnk_C = find( gridWindSc == windThresh );

%% export

intParam.intElem = intElem;
intParam.intConn = intConn;
intParam.cutLen = cutLen;
intParam.elemNorms = elemNorms;
intParam.fracElem = fracElem;
intParam.fracFaces = fracFaces;
intParam.intNodes = intNodes;
intParam.intPoints = intPoints;
intParam.intNorms = intNorms;
intParam.elTypes = elTypes;
intParam.exclNodes_U = exclNodes_U;
intParam.exclNodes_T = exclNodes_T;
intParam.intNodesUnk_U = intNodesUnk_U;
intParam.intNodesUnk_T = intNodesUnk_T;
intParam.intNodesUnk_C = intNodesUnk_C;

end

