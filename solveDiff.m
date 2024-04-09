function [ conc ] = solveDiff( ceq, matPar, geomPar, numPar, intParam )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

alpD = matPar.alpD;
ksD = matPar.ksD;

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dh = geomPar.dh;

BCdiff = geomPar.BCdiff;

intElem = intParam.intElem;
intConn = intParam.intConn;
cutLen = intParam.cutLen;
fracElem = intParam.fracElem;
fracFaces = intParam.fracFaces;
intNodes = intParam.intNodes;
elTypes = intParam.elTypes;
exclNodes_T = intParam.exclNodes_T;

KR = zeros( Nx*Ny, 36 );
KC = zeros( Nx*Ny, 36 );
KV = zeros( Nx*Ny, 36 );
b = zeros( Nx*Ny, 1 );

for ii=1:((Nx-1)*(Ny-1))

    ind_bl = ii + floor((ii-0.5)/(Nx-1));
    ind_br = ii + 1 + floor((ii-0.5)/(Nx-1));
    ind_tl = ii + Nx + floor((ii-0.5)/(Nx-1));
    ind_tr = ii + Nx + 1 + floor((ii-0.5)/(Nx-1));
    
    elem = [ ind_tl  ind_bl  ind_br   1 ;
             ind_br  ind_tr  ind_tl  -1 ];
         
    elem_inds = [ 2*ii-1 ;
                  2*ii   ];

    %. bottom el.
    gradva{1} = [  0  1 ;
                  -1 -1 ;
                   1  0 ];
    %. top el.
    gradva{2} = [  0 -1 ;
                   1  1 ;
                  -1  0 ];
              
    for mm=1:2
        
        elType_loc = elTypes( elem_inds(mm) );
        if ( elType_loc == 2 )
            %. untransformed
            area = 0;
        elseif ( elType_loc == 1 )
            %. transformed
            area = 1;
        else
            %. intersected
            indEl_loc = find( intElem == elem_inds(mm), 1, 'first' );
            area = 1 - fracElem(indEl_loc);
        end

        gradv = gradva{mm};
        for jj=1:3
            for kk=1:3
                KR( ii, 18*(mm-1)+3*(jj-1)+kk ) = elem(mm,jj);
                KC( ii, 18*(mm-1)+3*(jj-1)+kk ) = elem(mm,kk);
                KV( ii, 18*(mm-1)+3*(jj-1)+kk ) = area*gradv(jj,:)*gradv(kk,:).'*(1/2);
            end
        end
        
        %. line integrals
        if ( elType_loc == 3 )
            
            %. intersected points
            indEl_loc = find( intElem == elem_inds(mm), 1, 'first' );
            elIntPts = intConn(indEl_loc,:);
            elInt_p1n1ind = intNodes( elIntPts(1), 1 );
            elInt_p1n1w = fracFaces( elIntPts(1) );
            elInt_p1n2ind = intNodes( elIntPts(1), 2 );
            elInt_p1n2w = 1 - fracFaces( elIntPts(1) );
            elInt_p2n1ind = intNodes( elIntPts(2), 1 );
            elInt_p2n1w = fracFaces( elIntPts(2) );
            elInt_p2n2ind = intNodes( elIntPts(2), 2 );
            elInt_p2n2w = 1 - fracFaces( elIntPts(2) );

            for jj=1:3
                for kk=1:3
                    node1 = elem(mm,jj);
                    node2 = elem(mm,kk);
                    
                    phi1_p1n1 = 0;
                    if ( elInt_p1n1ind == node1 )
                        phi1_p1n1 = 1;
                    end
                    phi1_p1n2 = 0;
                    if ( elInt_p1n2ind == node1 )
                        phi1_p1n2 = 1;
                    end
                    phi1_p2n1 = 0;
                    if ( elInt_p2n1ind == node1 )
                        phi1_p2n1 = 1;
                    end
                    phi1_p2n2 = 0;
                    if ( elInt_p2n2ind == node1 )
                        phi1_p2n2 = 1;
                    end
                    phi1_p1 = phi1_p1n1 * elInt_p1n1w + phi1_p1n2 * elInt_p1n2w;
                    phi1_p2 = phi1_p2n1 * elInt_p2n1w + phi1_p2n2 * elInt_p2n2w;
                    
                    phi2_p1n1 = 0;
                    if ( elInt_p1n1ind == node2 )
                        phi2_p1n1 = 1;
                    end
                    phi2_p1n2 = 0;
                    if ( elInt_p1n2ind == node2 )
                        phi2_p1n2 = 1;
                    end
                    phi2_p2n1 = 0;
                    if ( elInt_p2n1ind == node2 )
                        phi2_p2n1 = 1;
                    end
                    phi2_p2n2 = 0;
                    if ( elInt_p2n2ind == node2 )
                        phi2_p2n2 = 1;
                    end
                    phi2_p1 = phi2_p1n1 * elInt_p1n1w + phi2_p1n2 * elInt_p1n2w;
                    phi2_p2 = phi2_p2n1 * elInt_p2n1w + phi2_p2n2 * elInt_p2n2w;
                    
                    secLen = cutLen( indEl_loc );
                    
                    %. integral, LHS
                    int = ksD * secLen * ( 2*phi1_p1*phi2_p1 + phi1_p1*phi2_p2 + phi1_p2*phi2_p1 + 2*phi1_p2*phi2_p2 ) * (1/6);
                    KR( ii, 9+18*(mm-1)+3*(jj-1)+kk ) = node1;
                    KC( ii, 9+18*(mm-1)+3*(jj-1)+kk ) = node2;
                    KV( ii, 9+18*(mm-1)+3*(jj-1)+kk ) = int;
                    
                    %. integral, RHS
                    if (jj == kk)
                        ceq_p1 = ceq( elIntPts(1) );
                        ceq_p2 = ceq( elIntPts(2) );
                        int = ksD * secLen * ( 2*phi1_p1*ceq_p1 + phi1_p1*ceq_p2 + phi1_p2*ceq_p1 + 2*phi1_p2*ceq_p2 ) * (1/6);
                        
                        b( node1 ) = b( node1 ) + int;
                    end
                end
            end
        end
    end
end

%. stablisation integrals
nS = size(intElem,1)*3;
GR = zeros(nS,16);
GC = zeros(nS,16);
GV = zeros(nS,16);
indS = 1;
kap = numPar.STAB_PAR;
for ii=1:(Nx-1)
    for jj=1:(Ny-1)

        %. squares
        sq_ar = ii + (jj-1)*(Nx-1);
        sq_al = ii + (jj-1)*(Nx-1) - 1;
        sq_br = ii + (jj-2)*(Nx-1);

        for mm=1:3

            if (mm==1)&&(ii>=2)
                stabCond = 1;
            elseif (mm==2)&&(jj>=2)
                stabCond = 1;
            elseif (mm==3)
                stabCond = 1;
            else
                stabCond = 0;
            end

            if ( stabCond == 1 )

                if (mm==1)
                    %. vertical edge
                    el_L = 2*sq_al;
                    el_R = 2*sq_ar-1;
                    n1 = ii + jj*Nx - 1;
                    n2 = ii + jj*Nx;
                    n3 = ii + (jj-1)*Nx;
                    n4 = ii + (jj-1)*Nx + 1;
                    len = 1;
                elseif (mm==2)
                    %. horizontal edge
                    el_L = 2*sq_br;
                    el_R = 2*sq_ar-1;
                    n1 = ii + jj*Nx;
                    n2 = ii + (jj-1)*Nx + 1;
                    n3 = ii + (jj-1)*Nx;
                    n4 = ii + (jj-2)*Nx + 1;
                    len = 1;
                elseif (mm==3)
                    %. diagonal edge
                    el_L = 2*sq_ar;
                    el_R = 2*sq_ar-1;
                    n1 = ii + (jj-1)*Nx;
                    n2 = ii + jj*Nx;
                    n3 = ii + (jj-1)*Nx + 1;
                    n4 = ii + jj*Nx + 1;
                    len = 2*sqrt(2);
                end
        
                if ( elTypes(el_L) == 3 )||( elTypes(el_R) == 3 )
                    %. stabilise edge
                    ir = [ n1 n1 n1 n1 n2 n2 n2 n2 n3 n3 n3 n3 n4 n4 n4 n4 ];
                    ic = [ n1 n2 n3 n4 n1 n2 n3 n4 n1 n2 n3 n4 n1 n2 n3 n4 ];
                    vs = [  1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1 ];

                    if ( elTypes(el_L) ~= 2 )&&( elTypes(el_R) ~= 2 )
                        %. transformed or intersected
                        GR(indS,:) = ir;
                        GC(indS,:) = ic;
                        GV(indS,:) = vs*kap*len;
                        indS = indS + 1;
                    end
                end
            end
        end
    end
end

GRc = GR(1:(indS-1),:);
GCc = GC(1:(indS-1),:);
GVc = GV(1:(indS-1),:);
GRr = GRc(:);
GCr = GCc(:);
GVr = GVc(:);
Gs = sparse( GRr, GCr, GVr, Nx*Ny, Nx*Ny );

KRr = KR(:);
KCr = KC(:);
KVr = KV(:);
idxrem = find( KRr==0 );
KRr( idxrem ) = [];
KCr( idxrem ) = [];
KVr( idxrem ) = [];
Ks = sparse(KRr,KCr,KVr);
Ks = Ks + Gs;

%. prescribe unused nodes
Kexcl = sparse( exclNodes_T, exclNodes_T, 1, Nx*Ny, Nx*Ny );
Ks = Ks + Kexcl;

if ( BCdiff == 0 )||( BCdiff == 1 )
    %. bottom side - mixed B.C.
    KaddR = zeros(Nx-1,4);
    KaddC = zeros(Nx-1,4);
    KaddV = zeros(Nx-1,4);
    for ii=1:(Nx-1)

        ind_bl = ii + floor((ii-0.5)/(Nx-1));
        ind_br = ii + 1 + floor((ii-0.5)/(Nx-1));

        KaddR( ii, 1 ) = ind_bl;
        KaddC( ii, 1 ) = ind_bl;
        KaddV( ii, 1 ) = alpD*(1/3)*dh;
        KaddR( ii, 2 ) = ind_br;
        KaddC( ii, 2 ) = ind_br;
        KaddV( ii, 2 ) = alpD*(1/3)*dh;
        KaddR( ii, 3 ) = ind_bl;
        KaddC( ii, 3 ) = ind_br;
        KaddV( ii, 3 ) = alpD*(1/6)*dh;
        KaddR( ii, 4 ) = ind_br;
        KaddC( ii, 4 ) = ind_bl;
        KaddV( ii, 4 ) = alpD*(1/6)*dh;
        b( ind_bl ) = b( ind_bl ) + alpD*dh/2;
        b( ind_br ) = b( ind_br ) + alpD*dh/2;

    end
    Kadd = sparse( KaddR, KaddC, KaddV, Nx*Ny, Nx*Ny );
    Ks = Ks + Kadd;
end

if ( BCdiff == 1 )
    %. top side - mixed B.C.
    KaddR = zeros(Nx-1,4);
    KaddC = zeros(Nx-1,4);
    KaddV = zeros(Nx-1,4);
    for ii=1:(Nx-1)

        sqind = (Nx-1)*(Ny-2) + ii;

        ind_tl = sqind + Nx + floor((sqind-0.5)/(Nx-1));
        ind_tr = sqind + Nx + 1 + floor((sqind-0.5)/(Nx-1));

        KaddR( ii, 1 ) = ind_tl;
        KaddC( ii, 1 ) = ind_tl;
        KaddV( ii, 1 ) = alpD*(1/3)*dh;
        KaddR( ii, 2 ) = ind_tr;
        KaddC( ii, 2 ) = ind_tr;
        KaddV( ii, 2 ) = alpD*(1/3)*dh;
        KaddR( ii, 3 ) = ind_tl;
        KaddC( ii, 3 ) = ind_tr;
        KaddV( ii, 3 ) = alpD*(1/6)*dh;
        KaddR( ii, 4 ) = ind_tr;
        KaddC( ii, 4 ) = ind_tl;
        KaddV( ii, 4 ) = alpD*(1/6)*dh;
        b( ind_tl ) = b( ind_tl ) + alpD*dh/2;
        b( ind_tr ) = b( ind_tr ) + alpD*dh/2;

    end
    Kadd = sparse( KaddR, KaddC, KaddV, Nx*Ny, Nx*Ny );
    Ks = Ks + Kadd;

    %. left side - mixed B.C.
    KaddR = zeros(Ny-1,4);
    KaddC = zeros(Ny-1,4);
    KaddV = zeros(Ny-1,4);
    for ii=1:(Ny-1)

        sqind = 1 + (ii-1)*(Nx-1);

        ind_bl = sqind + floor((sqind-0.5)/(Nx-1));
        ind_tl = sqind + Nx + floor((sqind-0.5)/(Nx-1));

        KaddR( ii, 1 ) = ind_bl;
        KaddC( ii, 1 ) = ind_bl;
        KaddV( ii, 1 ) = alpD*(1/3)*dh;
        KaddR( ii, 2 ) = ind_tl;
        KaddC( ii, 2 ) = ind_tl;
        KaddV( ii, 2 ) = alpD*(1/3)*dh;
        KaddR( ii, 3 ) = ind_bl;
        KaddC( ii, 3 ) = ind_tl;
        KaddV( ii, 3 ) = alpD*(1/6)*dh;
        KaddR( ii, 4 ) = ind_tl;
        KaddC( ii, 4 ) = ind_bl;
        KaddV( ii, 4 ) = alpD*(1/6)*dh;
        b( ind_bl ) = b( ind_bl ) + alpD*dh/2;
        b( ind_tl ) = b( ind_tl ) + alpD*dh/2;

    end
    Kadd = sparse( KaddR, KaddC, KaddV, Nx*Ny, Nx*Ny );
    Ks = Ks + Kadd;

    %. right side - mixed B.C.
    KaddR = zeros(Ny-1,4);
    KaddC = zeros(Ny-1,4);
    KaddV = zeros(Ny-1,4);
    for ii=1:(Ny-1)

        sqind = (Nx-1) + (ii-1)*(Nx-1);

        ind_br = sqind + 1 + floor((sqind-0.5)/(Nx-1));
        ind_tr = sqind + Nx + 1 + floor((sqind-0.5)/(Nx-1));

        KaddR( ii, 1 ) = ind_tr;
        KaddC( ii, 1 ) = ind_tr;
        KaddV( ii, 1 ) = alpD*(1/3)*dh;
        KaddR( ii, 2 ) = ind_br;
        KaddC( ii, 2 ) = ind_br;
        KaddV( ii, 2 ) = alpD*(1/3)*dh;
        KaddR( ii, 3 ) = ind_tr;
        KaddC( ii, 3 ) = ind_br;
        KaddV( ii, 3 ) = alpD*(1/6)*dh;
        KaddR( ii, 4 ) = ind_br;
        KaddC( ii, 4 ) = ind_tr;
        KaddV( ii, 4 ) = alpD*(1/6)*dh;
        b( ind_tr ) = b( ind_tr ) + alpD*dh/2;
        b( ind_br ) = b( ind_br ) + alpD*dh/2;

    end
    Kadd = sparse( KaddR, KaddC, KaddV, Nx*Ny, Nx*Ny );
    Ks = Ks + Kadd;
end

%. solve
conc = Ks \ b;

end

