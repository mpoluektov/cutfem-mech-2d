function [ F, Js, resElem ] = calcMechResPhaseTran( u, matPar, geomPar, numPar, intParam, prevState, g, dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

elGU = matPar.elGU;
elGT = matPar.elGT;
elKU = matPar.elKU;
elKT = matPar.elKT;
etaU = matPar.etaU;
etaT = matPar.etaT;

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dh = geomPar.dh;

F_U_od = zeros( Nx*Ny, 1 );
F_U_ev = zeros( Nx*Ny, 1 );
F_T_od = zeros( Nx*Ny, 1 );
F_T_ev = zeros( Nx*Ny, 1 );

uU = u( 1:(2*Nx*Ny) );
uT = u( (2*Nx*Ny+1):(4*Nx*Ny) );

JR = zeros( Nx*Ny, 36 );
JC = zeros( Nx*Ny, 36 );
JVx_U = zeros( Nx*Ny, 36 );
JVy_U = zeros( Nx*Ny, 36 );
JVx_T = zeros( Nx*Ny, 36 );
JVy_T = zeros( Nx*Ny, 36 );
JVx_U_ex = zeros( Nx*Ny, 36 );
JVy_U_ex = zeros( Nx*Ny, 36 );
JVx_T_ex = zeros( Nx*Ny, 36 );
JVy_T_ex = zeros( Nx*Ny, 36 );
JR_n = zeros( Nx*Ny, 72 );
JC_n = zeros( Nx*Ny, 72 );
JVx_U_n = zeros( Nx*Ny, 72 );
JVy_U_n = zeros( Nx*Ny, 72 );
JVx_T_n = zeros( Nx*Ny, 72 );
JVy_T_n = zeros( Nx*Ny, 72 );

%% get interface elements and nodes

intElem = intParam.intElem;
intConn = intParam.intConn;
cutLen = intParam.cutLen;
elemNorms = intParam.elemNorms;
fracElem = intParam.fracElem;
fracFaces = intParam.fracFaces;
intNodes = intParam.intNodes;
elTypes = intParam.elTypes;
exclNodes_U = intParam.exclNodes_U;
exclNodes_T = intParam.exclNodes_T;
intNodesUnk_C = intParam.intNodesUnk_C;

%% previous plastic def.

Cp_prev_U_all = prevState.Cp_prev_U_all;
Cp_prev_T_all = prevState.Cp_prev_T_all;
Cp_prev_U_all_ev = Cp_prev_U_all( 2:2:(2*(Nx-1)*(Ny-1)), : );
Cp_prev_U_all_od = Cp_prev_U_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : );
Cp_prev_T_all_ev = Cp_prev_T_all( 2:2:(2*(Nx-1)*(Ny-1)), : );
Cp_prev_T_all_od = Cp_prev_T_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : );

%% calculate Piola-Kirchhoff stresses of elements

Pe_U_all_od = zeros( (Nx-1)*(Ny-1), 9 );
Pe_U_all_ev = zeros( (Nx-1)*(Ny-1), 9 );
d_PeT_d_Fe_U_all_od = zeros( (Nx-1)*(Ny-1), 81 );
d_PeT_d_Fe_U_all_ev = zeros( (Nx-1)*(Ny-1), 81 );
Fe_U_all_od = zeros( (Nx-1)*(Ny-1), 9 );
Fe_U_all_ev = zeros( (Nx-1)*(Ny-1), 9 );
We_U_all_od = zeros( (Nx-1)*(Ny-1), 1 );
We_U_all_ev = zeros( (Nx-1)*(Ny-1), 1 );
Cp_U_all_od = zeros( (Nx-1)*(Ny-1), 9 );
Cp_U_all_ev = zeros( (Nx-1)*(Ny-1), 9 );
Pe_T_all_od = zeros( (Nx-1)*(Ny-1), 9 );
Pe_T_all_ev = zeros( (Nx-1)*(Ny-1), 9 );
d_PeT_d_Fe_T_all_od = zeros( (Nx-1)*(Ny-1), 81 );
d_PeT_d_Fe_T_all_ev = zeros( (Nx-1)*(Ny-1), 81 );
Fe_T_all_od = zeros( (Nx-1)*(Ny-1), 9 );
Fe_T_all_ev = zeros( (Nx-1)*(Ny-1), 9 );
We_T_all_od = zeros( (Nx-1)*(Ny-1), 1 );
We_T_all_ev = zeros( (Nx-1)*(Ny-1), 1 );
Cp_T_all_od = zeros( (Nx-1)*(Ny-1), 9 );
Cp_T_all_ev = zeros( (Nx-1)*(Ny-1), 9 );
elemNodeInds_all_od = zeros( (Nx-1)*(Ny-1), 4 );
elemNodeInds_all_ev = zeros( (Nx-1)*(Ny-1), 4 );

% parfor ii=1:((Nx-1)*(Ny-1))
for ii=1:((Nx-1)*(Ny-1))

    ind_bl = ii + floor((ii-0.5)/(Nx-1));
    ind_br = ii + 1 + floor((ii-0.5)/(Nx-1));
    ind_tl = ii + Nx + floor((ii-0.5)/(Nx-1));
    ind_tr = ii + Nx + 1 + floor((ii-0.5)/(Nx-1));
    
    elem = [ ind_tl  ind_bl  ind_br   1 ;
             ind_br  ind_tr  ind_tl  -1 ];
    
    for jj=1:2
        
        ind1x = 2*elem(jj,1)-1;
        ind1y = 2*elem(jj,1);
        ind2x = 2*elem(jj,2)-1;
        ind2y = 2*elem(jj,2);
        ind3x = 2*elem(jj,3)-1;
        ind3y = 2*elem(jj,3);
        
        elInd = 2*ii+jj-2;
        
        if ( elTypes(elInd) == 2 )||( elTypes(elInd) == 3 )
            %. untransformed or intersected
            ud = elem(jj,4) * [ uU(ind3x)-uU(ind2x)  uU(ind1x)-uU(ind2x) ;
                                uU(ind3y)-uU(ind2y)  uU(ind1y)-uU(ind2y) ] * (1/dh);
            ude = zeros(3);
            ude(1:2,1:2) = ud;
            Fe = ude + eye(3);

            if ( jj == 1 )
                Cp_prev = vec2tens(Cp_prev_U_all_od(ii,:));
            else
                Cp_prev = vec2tens(Cp_prev_U_all_ev(ii,:));
            end

            [ Pe, d_PeT_d_Fe, We, Cp ] = constLaw( Fe, Cp_prev, elKU, elGU, etaU, 1, dt );

            if ( jj == 1 )
                Pe_U_all_od(ii,:) = tens2vec(Pe).';
                d_PeT_d_Fe_U_all_od(ii,:) = d_PeT_d_Fe(:);
                Fe_U_all_od(ii,:) = tens2vec(Fe).';
                We_U_all_od(ii,:) = We;
                Cp_U_all_od(ii,:) = tens2vec(Cp).';
            else
                Pe_U_all_ev(ii,:) = tens2vec(Pe).';
                d_PeT_d_Fe_U_all_ev(ii,:) = d_PeT_d_Fe(:);
                Fe_U_all_ev(ii,:) = tens2vec(Fe).';
                We_U_all_ev(ii,:) = We;
                Cp_U_all_ev(ii,:) = tens2vec(Cp).';
            end
        end
        if ( elTypes(elInd) == 1 )||( elTypes(elInd) == 3 )
            %. transformed or intersected
            ud = elem(jj,4) * [ uT(ind3x)-uT(ind2x)  uT(ind1x)-uT(ind2x) ;
                                uT(ind3y)-uT(ind2y)  uT(ind1y)-uT(ind2y) ] * (1/dh);
            ude = zeros(3);
            ude(1:2,1:2) = ud;
            Fe = ude + eye(3);

            if ( jj == 1 )
                Cp_prev = vec2tens(Cp_prev_T_all_od(ii,:));
            else
                Cp_prev = vec2tens(Cp_prev_T_all_ev(ii,:));
            end

            [ Pe, d_PeT_d_Fe, We, Cp ] = constLaw( Fe, Cp_prev, elKT, elGT, etaT, g, dt );

            if ( jj == 1 )
                Pe_T_all_od(ii,:) = tens2vec(Pe).';
                d_PeT_d_Fe_T_all_od(ii,:) = d_PeT_d_Fe(:);
                Fe_T_all_od(ii,:) = tens2vec(Fe).';
                We_T_all_od(ii,:) = We;
                Cp_T_all_od(ii,:) = tens2vec(Cp).';
            else
                Pe_T_all_ev(ii,:) = tens2vec(Pe).';
                d_PeT_d_Fe_T_all_ev(ii,:) = d_PeT_d_Fe(:);
                Fe_T_all_ev(ii,:) = tens2vec(Fe).';
                We_T_all_ev(ii,:) = We;
                Cp_T_all_ev(ii,:) = tens2vec(Cp).';
            end
        end
        if ( jj == 1 )
            elemNodeInds_all_od(ii,:) = elem(jj,:);
        else
            elemNodeInds_all_ev(ii,:) = elem(jj,:);
        end
    end
end

Pe_U_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = Pe_U_all_ev;
Pe_U_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = Pe_U_all_od;
d_PeT_d_Fe_U_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d_PeT_d_Fe_U_all_ev;
d_PeT_d_Fe_U_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d_PeT_d_Fe_U_all_od;
Fe_U_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = Fe_U_all_ev;
Fe_U_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = Fe_U_all_od;
We_U_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = We_U_all_ev;
We_U_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = We_U_all_od;
Cp_U_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = Cp_U_all_ev;
Cp_U_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = Cp_U_all_od;
Pe_T_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = Pe_T_all_ev;
Pe_T_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = Pe_T_all_od;
d_PeT_d_Fe_T_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = d_PeT_d_Fe_T_all_ev;
d_PeT_d_Fe_T_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = d_PeT_d_Fe_T_all_od;
Fe_T_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = Fe_T_all_ev;
Fe_T_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = Fe_T_all_od;
We_T_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = We_T_all_ev;
We_T_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = We_T_all_od;
Cp_T_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = Cp_T_all_ev;
Cp_T_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = Cp_T_all_od;
elemNodeInds_all( 2:2:(2*(Nx-1)*(Ny-1)), : ) = elemNodeInds_all_ev;
elemNodeInds_all( 1:2:(2*(Nx-1)*(Ny-1)-1), : ) = elemNodeInds_all_od;

%% export

resElem.Pe_U_all = Pe_U_all;
resElem.Pe_T_all = Pe_T_all;
resElem.d_PeT_d_Fe_U_all = d_PeT_d_Fe_U_all;
resElem.d_PeT_d_Fe_T_all = d_PeT_d_Fe_T_all;
resElem.Fe_U_all = Fe_U_all;
resElem.Fe_T_all = Fe_T_all;
resElem.We_U_all = We_U_all;
resElem.We_T_all = We_T_all;
resElem.Cp_U_all = Cp_U_all;
resElem.Cp_T_all = Cp_T_all;

%% calculate integrals

gradv = [  1   0 ;
          -1  -1 ;
           0   1 ;
           0  -1 ;
          -1   0 ;
           1   1 ;
           0   0 ; 
           0   0 ; 
           0   0 ; 
           0   0 ; 
           0   0 ; 
           0   0 ];
d_Fe_d_u_r = [ 0  0 -1  0  1  0 ;
               0  1  0 -1  0  0 ;
               0  0  0  0  0  0 ;
               1  0 -1  0  0  0 ;
               0  0  0  0  0  0 ;
               0  0  0  0  0  0 ;
               0  0  0 -1  0  1 ;
               0  0  0  0  0  0 ;
               0  0  0  0  0  0 ];
elem_edeges = [  1  7  1 ;
                 1  6  2 ;
                 4  1  3 ;
                 8  4  2 ;
                 2  4  1 ;
                 9  2  3 ;
                 2  5  2 ;
                10  5  1 ;
                 5  3  3 ;
                 3 11  2 ;
                 3  6  1 ;
                 6 12  3 ];

% parfor ii=1:(Nx*Ny)
for ii=1:(Nx*Ny)
               
    %. squares
    sq_ar = ii - floor( (ii-0.5)/Nx );
    sq_al = ii - floor( (ii-0.5)/Nx ) - 1;
    sq_br = ii - floor( (ii-0.5)/Nx ) - (Nx-1);
    sq_bl = ii - floor( (ii-0.5)/Nx ) - (Nx-1) - 1;

    sq_all = ii - floor( (ii-0.5)/Nx ) - 2;
    sq_aal = ii - floor( (ii-0.5)/Nx ) + (Nx-1) - 1;
    sq_brr = ii - floor( (ii-0.5)/Nx ) - (Nx-1) + 1;
    sq_bbr = ii - floor( (ii-0.5)/Nx ) - 2*(Nx-1);
    
    elem_inds = [ 2*sq_al-1 ;
                  2*sq_ar-1 ;
                  2*sq_br-1 ;
                    2*sq_al ;
                    2*sq_br ;
                    2*sq_bl ;
                   2*sq_all ;
                 2*sq_aal-1 ;
                    2*sq_ar ;
                 2*sq_brr-1 ;
                   2*sq_bbr ;
                  2*sq_bl-1 ];

    if ( ii == 1 )
        %. LB corner
        inda = 2;
    elseif ( ii == Nx*Ny-Nx+1 )
        %. LT corner
        inda = [ 3 5 ];
    elseif ( ii == Nx )
        %. RB corner
        inda = [ 1 4 ];
    elseif ( ii == Nx*Ny )
        %. RT corner
        inda = 6;
    elseif ( rem(ii,Nx) == 1 )
        %. L side
        inda = [ 2 3 5 ];
    elseif ( ii > Nx*Ny-Nx+1 )
        %. T side
        inda = [ 3 5 6 ];
    elseif ( rem(ii,Nx) == 0 )
        %. R side
        inda = [ 1 4 6 ];
    elseif ( ii < Nx )
        %. B side
        inda = [ 1 2 4 ];
    else
        %. inside
        inda = [ 1 2 3 4 5 6 ];
    end
    
    %. find edges
    inde = 1:12;
    if ( rem(ii,Nx) == 1 )
        %. L side
        inde( inde==1 ) = [];
        inde( inde==2 ) = [];
        inde( inde==3 ) = [];
        inde( inde==4 ) = [];
        inde( inde==5 ) = [];
        inde( inde==11 ) = [];
        inde( inde==12 ) = [];
    end
    if ( rem(ii,Nx) == 2 )
        %. L side dist. by h
        inde( inde==1 ) = [];
    end
    if ( rem(ii,Nx) == 0 )
        %. R side
        inde( inde==5 ) = [];
        inde( inde==6 ) = [];
        inde( inde==7 ) = [];
        inde( inde==8 ) = [];
        inde( inde==9 ) = [];
        inde( inde==10 ) = [];
        inde( inde==11 ) = [];
    end
    if ( rem(ii,Nx) == Nx-1 )
        %. R side dist. by h
        inde( inde==8 ) = [];
    end
    if ( ii <= Nx )
        %. B side
        inde( inde==2 ) = [];
        inde( inde==7 ) = [];
        inde( inde==8 ) = [];
        inde( inde==9 ) = [];
        inde( inde==10 ) = [];
        inde( inde==11 ) = [];
        inde( inde==12 ) = [];
    end
    if ( ii <= 2*Nx )&&( ii >= Nx+1 )
        %. B side dist. by h
        inde( inde==10 ) = [];
    end
    if ( ii >= Nx*Ny-Nx+1 )
        %. T side
        inde( inde==1 ) = [];
        inde( inde==2 ) = [];
        inde( inde==3 ) = [];
        inde( inde==4 ) = [];
        inde( inde==5 ) = [];
        inde( inde==6 ) = [];
        inde( inde==7 ) = [];
    end
    if ( ii >= Nx*Ny-2*Nx+1 )&&( ii <= Nx*Ny-Nx )
        %. T side dist. by h
        inde( inde==4 ) = [];
    end
    
    TR = zeros( 1, 36 );
    TC = zeros( 1, 36 );
    TVx_U = zeros( 1, 36 );
    TVy_U = zeros( 1, 36 );
    TVx_T = zeros( 1, 36 );
    TVy_T = zeros( 1, 36 );
    TVx_U_ex = zeros( 1, 36 );
    TVy_U_ex = zeros( 1, 36 );
    TVx_T_ex = zeros( 1, 36 );
    TVy_T_ex = zeros( 1, 36 );
    TR_n = zeros( 1, 72 );
    TC_n = zeros( 1, 72 );
    TVx_U_n = zeros( 1, 72 );
    TVy_U_n = zeros( 1, 72 );
    TVx_T_n = zeros( 1, 72 );
    TVy_T_n = zeros( 1, 72 );

    for jj=1:size(inda,2)
        
        Pe_U = vec2tens( Pe_U_all( elem_inds(inda(jj)), : ) );
        Pe_T = vec2tens( Pe_T_all( elem_inds(inda(jj)), : ) );
        Pec_U = Pe_U(1:2,1:2);
        Pec_T = Pe_T(1:2,1:2);
        
        elemNodeInds = elemNodeInds_all( elem_inds(inda(jj)), : );
        
        elType_loc = elTypes( elem_inds(inda(jj)) );
        if ( elType_loc == 2 )
            %. untransformed
            area_U = 1;
            area_T = 0;
        elseif ( elType_loc == 1 )
            %. transformed
            area_U = 0;
            area_T = 1;
        else
            %. intersected
            indEl_loc = find( intElem == elem_inds(inda(jj)), 1, 'first' );
            area_U = fracElem(indEl_loc);
            area_T = 1 - fracElem(indEl_loc);
        end
        
        gphi1 = (gradv(inda(jj),:).') * (1/dh) * [1 0];
        gphi2 = (gradv(inda(jj),:).') * (1/dh) * [0 1];
        
        int1_U = area_U * tensTrace2( Pec_U * gphi1 ) * (1/2) * (dh^2);
        int2_U = area_U * tensTrace2( Pec_U * gphi2 ) * (1/2) * (dh^2);
        
        F_U_od(ii) = F_U_od(ii) + int1_U;
        F_U_ev(ii) = F_U_ev(ii) + int2_U;
        
        int1_T = area_T * tensTrace2( Pec_T * gphi1 ) * (1/2) * (dh^2);
        int2_T = area_T * tensTrace2( Pec_T * gphi2 ) * (1/2) * (dh^2);
        
        F_T_od(ii) = F_T_od(ii) + int1_T;
        F_T_ev(ii) = F_T_ev(ii) + int2_T;
        
        %. N.R. derivatives
        d_PeT_d_Fe_U_r = d_PeT_d_Fe_U_all( elem_inds(inda(jj)), : );
        d_PeT_d_Fe_U = reshape( d_PeT_d_Fe_U_r, 9, 9 );
        d_PeT_d_Fe_T_r = d_PeT_d_Fe_T_all( elem_inds(inda(jj)), : );
        d_PeT_d_Fe_T = reshape( d_PeT_d_Fe_T_r, 9, 9 );

        gphi1e = zeros(3);
        gphi1e(1:2,1:2) = gphi1;
        vgphi1e = tens2vec(gphi1e);
        gphi2e = zeros(3);
        gphi2e(1:2,1:2) = gphi2;
        vgphi2e = tens2vec(gphi2e);

        d_Fe_d_u = elemNodeInds(4) * d_Fe_d_u_r * (1/dh);
        
        d_int1_U_d_Fe = area_U * vgphi1e.' * d_PeT_d_Fe_U * (1/2) * (dh^2);
        d_int2_U_d_Fe = area_U * vgphi2e.' * d_PeT_d_Fe_U * (1/2) * (dh^2);
        d_int1_U_d_u = d_int1_U_d_Fe * d_Fe_d_u;
        d_int2_U_d_u = d_int2_U_d_Fe * d_Fe_d_u;

        d_int1_T_d_Fe = area_T * vgphi1e.' * d_PeT_d_Fe_T * (1/2) * (dh^2);
        d_int2_T_d_Fe = area_T * vgphi2e.' * d_PeT_d_Fe_T * (1/2) * (dh^2);
        d_int1_T_d_u = d_int1_T_d_Fe * d_Fe_d_u;
        d_int2_T_d_u = d_int2_T_d_Fe * d_Fe_d_u;

        ind1x = 2*elemNodeInds(1)-1;
        ind1y = 2*elemNodeInds(1);
        ind2x = 2*elemNodeInds(2)-1;
        ind2y = 2*elemNodeInds(2);
        ind3x = 2*elemNodeInds(3)-1;
        ind3y = 2*elemNodeInds(3);
        idx = [ ind1x ind1y ind2x ind2y ind3x ind3y ];
        TR( (6*jj-5):(6*jj) ) = repmat( ii, 1, 6 );
        TC( (6*jj-5):(6*jj) ) = idx;
        TVx_U( (6*jj-5):(6*jj) ) = d_int1_U_d_u;
        TVy_U( (6*jj-5):(6*jj) ) = d_int2_U_d_u;
        TVx_T( (6*jj-5):(6*jj) ) = d_int1_T_d_u;
        TVy_T( (6*jj-5):(6*jj) ) = d_int2_T_d_u;
        
        %. line integrals
        if ( elType_loc == 3 )
            
            %. intersected
            indEl_loc = find( intElem == elem_inds(inda(jj)), 1, 'first' );
            elIntPts = intConn(indEl_loc,:);
            elInt_p1n1ind = intNodes( elIntPts(1), 1 );
            elInt_p1n1w = fracFaces( elIntPts(1) );
            elInt_p1n2ind = intNodes( elIntPts(1), 2 );
            elInt_p1n2w = 1 - fracFaces( elIntPts(1) );
            elInt_p2n1ind = intNodes( elIntPts(2), 1 );
            elInt_p2n1w = fracFaces( elIntPts(2) );
            elInt_p2n2ind = intNodes( elIntPts(2), 2 );
            elInt_p2n2w = 1 - fracFaces( elIntPts(2) );
            phi_p1n1 = 0;
            if ( elInt_p1n1ind == ii )
                phi_p1n1 = 1;
            end
            phi_p1n2 = 0;
            if ( elInt_p1n2ind == ii )
                phi_p1n2 = 1;
            end
            phi_p2n1 = 0;
            if ( elInt_p2n1ind == ii )
                phi_p2n1 = 1;
            end
            phi_p2n2 = 0;
            if ( elInt_p2n2ind == ii )
                phi_p2n2 = 1;
            end
            phi_p1 = phi_p1n1 * elInt_p1n1w + phi_p1n2 * elInt_p1n2w;
            phi_p2 = phi_p2n1 * elInt_p2n1w + phi_p2n2 * elInt_p2n2w;
            
            uU_p1n1 = uU( (2*elInt_p1n1ind-1):(2*elInt_p1n1ind) );
            d_uU_p1n1_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p1n1ind-1) );
            d_uU_p1n1_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p1n1ind) );
            d_uU_p1n1_d_u(2,idx_loc) = 1;
            uU_p1n2 = uU( (2*elInt_p1n2ind-1):(2*elInt_p1n2ind) );
            d_uU_p1n2_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p1n2ind-1) );
            d_uU_p1n2_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p1n2ind) );
            d_uU_p1n2_d_u(2,idx_loc) = 1;
            uU_p2n1 = uU( (2*elInt_p2n1ind-1):(2*elInt_p2n1ind) );
            d_uU_p2n1_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p2n1ind-1) );
            d_uU_p2n1_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p2n1ind) );
            d_uU_p2n1_d_u(2,idx_loc) = 1;
            uU_p2n2 = uU( (2*elInt_p2n2ind-1):(2*elInt_p2n2ind) );
            d_uU_p2n2_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p2n2ind-1) );
            d_uU_p2n2_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p2n2ind) );
            d_uU_p2n2_d_u(2,idx_loc) = 1;
            uU_p1 = uU_p1n1 * elInt_p1n1w + uU_p1n2 * elInt_p1n2w;
            uU_p2 = uU_p2n1 * elInt_p2n1w + uU_p2n2 * elInt_p2n2w;
            d_uU_p1_d_u = d_uU_p1n1_d_u * elInt_p1n1w + d_uU_p1n2_d_u * elInt_p1n2w;
            d_uU_p2_d_u = d_uU_p2n1_d_u * elInt_p2n1w + d_uU_p2n2_d_u * elInt_p2n2w;
            
            uT_p1n1 = uT( (2*elInt_p1n1ind-1):(2*elInt_p1n1ind) );
            d_uT_p1n1_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p1n1ind-1) );
            d_uT_p1n1_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p1n1ind) );
            d_uT_p1n1_d_u(2,idx_loc) = 1;
            uT_p1n2 = uT( (2*elInt_p1n2ind-1):(2*elInt_p1n2ind) );
            d_uT_p1n2_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p1n2ind-1) );
            d_uT_p1n2_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p1n2ind) );
            d_uT_p1n2_d_u(2,idx_loc) = 1;
            uT_p2n1 = uT( (2*elInt_p2n1ind-1):(2*elInt_p2n1ind) );
            d_uT_p2n1_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p2n1ind-1) );
            d_uT_p2n1_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p2n1ind) );
            d_uT_p2n1_d_u(2,idx_loc) = 1;
            uT_p2n2 = uT( (2*elInt_p2n2ind-1):(2*elInt_p2n2ind) );
            d_uT_p2n2_d_u = zeros( 2, 6 );
            idx_loc = ( idx == (2*elInt_p2n2ind-1) );
            d_uT_p2n2_d_u(1,idx_loc) = 1;
            idx_loc = ( idx == (2*elInt_p2n2ind) );
            d_uT_p2n2_d_u(2,idx_loc) = 1;
            uT_p1 = uT_p1n1 * elInt_p1n1w + uT_p1n2 * elInt_p1n2w;
            uT_p2 = uT_p2n1 * elInt_p2n1w + uT_p2n2 * elInt_p2n2w;
            d_uT_p1_d_u = d_uT_p1n1_d_u * elInt_p1n1w + d_uT_p1n2_d_u * elInt_p1n2w;
            d_uT_p2_d_u = d_uT_p2n1_d_u * elInt_p2n1w + d_uT_p2n2_d_u * elInt_p2n2w;
            
            uJ_p1 = ( uT_p1 - uU_p1 );
            uJ_p2 = ( uT_p2 - uU_p2 );
            
            secLen = cutLen( indEl_loc );
            secNorm = elemNorms( indEl_loc, : ).';
            
            %. first integral: [phi].P.n
            int1 = secLen * (phi_p1 + phi_p2)*(1/2) * ( [1 0] * (Pec_U+Pec_T)*(1/2) * secNorm );
            int2 = secLen * (phi_p1 + phi_p2)*(1/2) * ( [0 1] * (Pec_U+Pec_T)*(1/2) * secNorm );
            %. untransformed
            F_U_od(ii) = F_U_od(ii) + int1;
            F_U_ev(ii) = F_U_ev(ii) + int2;
            %. transformed
            F_T_od(ii) = F_T_od(ii) - int1;
            F_T_ev(ii) = F_T_ev(ii) - int2;

            %. first integral, deriv.
            Nphi1e = zeros(3);
            Nphi1e(1:2,1:2) = secLen * (phi_p1 + phi_p2)*(1/2) * secNorm * [1 0] * (1/2);
            vNphi1e = tens2vec(Nphi1e);
            Nphi2e = zeros(3);
            Nphi2e(1:2,1:2) = secLen * (phi_p1 + phi_p2)*(1/2) * secNorm * [0 1] * (1/2);
            vNphi2e = tens2vec(Nphi2e);
            d_int1_d_Fe_U = vNphi1e.' * d_PeT_d_Fe_U;
            d_int2_d_Fe_U = vNphi2e.' * d_PeT_d_Fe_U;
            d_int1_d_uU = d_int1_d_Fe_U * d_Fe_d_u;
            d_int2_d_uU = d_int2_d_Fe_U * d_Fe_d_u;
            d_int1_d_Fe_T = vNphi1e.' * d_PeT_d_Fe_T;
            d_int2_d_Fe_T = vNphi2e.' * d_PeT_d_Fe_T;
            d_int1_d_uT = d_int1_d_Fe_T * d_Fe_d_u;
            d_int2_d_uT = d_int2_d_Fe_T * d_Fe_d_u;
            TVx_U( (6*jj-5):(6*jj) ) = TVx_U( (6*jj-5):(6*jj) ) + d_int1_d_uU;
            TVy_U( (6*jj-5):(6*jj) ) = TVy_U( (6*jj-5):(6*jj) ) + d_int2_d_uU;
            TVx_T( (6*jj-5):(6*jj) ) = TVx_T( (6*jj-5):(6*jj) ) - d_int1_d_uT;
            TVy_T( (6*jj-5):(6*jj) ) = TVy_T( (6*jj-5):(6*jj) ) - d_int2_d_uT;
            TVx_U_ex( (6*jj-5):(6*jj) ) = TVx_U_ex( (6*jj-5):(6*jj) ) + d_int1_d_uT;
            TVy_U_ex( (6*jj-5):(6*jj) ) = TVy_U_ex( (6*jj-5):(6*jj) ) + d_int2_d_uT;
            TVx_T_ex( (6*jj-5):(6*jj) ) = TVx_T_ex( (6*jj-5):(6*jj) ) - d_int1_d_uU;
            TVy_T_ex( (6*jj-5):(6*jj) ) = TVy_T_ex( (6*jj-5):(6*jj) ) - d_int2_d_uU;
            
            %. second integral: n.grad_phi.[u]
            %. gradients
            vgphi1T = tens2vec(gphi1e.');
            vgphi2T = tens2vec(gphi2e.');
            %. derivatives
            derP_U_gphi1 = vec2tens( d_PeT_d_Fe_U * vgphi1T );
            derP_U_gphi2 = vec2tens( d_PeT_d_Fe_U * vgphi2T );
            derP_T_gphi1 = vec2tens( d_PeT_d_Fe_T * vgphi1T );
            derP_T_gphi2 = vec2tens( d_PeT_d_Fe_T * vgphi2T );
            %. integrals
            int1_U = secLen * (secNorm.') * derP_U_gphi1(1:2,1:2) * (1/2) * (uJ_p1 + uJ_p2)*(1/2);
            int2_U = secLen * (secNorm.') * derP_U_gphi2(1:2,1:2) * (1/2) * (uJ_p1 + uJ_p2)*(1/2);
            int1_T = secLen * (secNorm.') * derP_T_gphi1(1:2,1:2) * (1/2) * (uJ_p1 + uJ_p2)*(1/2);
            int2_T = secLen * (secNorm.') * derP_T_gphi2(1:2,1:2) * (1/2) * (uJ_p1 + uJ_p2)*(1/2);
            %. untransformed
            F_U_od(ii) = F_U_od(ii) - int1_U;
            F_U_ev(ii) = F_U_ev(ii) - int2_U;
            %. transformed
            F_T_od(ii) = F_T_od(ii) - int1_T;
            F_T_ev(ii) = F_T_ev(ii) - int2_T;
            
            %. second integral, deriv.
            Fe_U = vec2tens( Fe_U_all(elem_inds(inda(jj)),:) );
            Fe_T = vec2tens( Fe_T_all(elem_inds(inda(jj)),:) );
            Cp_U = vec2tens( Cp_U_all(elem_inds(inda(jj)),:) );
            Cp_T = vec2tens( Cp_T_all(elem_inds(inda(jj)),:) );
            Am = secLen * (1/2) * (uJ_p1 + uJ_p2)*(1/2) * secNorm.';
            Ame = zeros(3);
            Ame(1:2,1:2) = Am;
            d_int1_U = secLen * (secNorm.') * derP_U_gphi1(1:2,1:2) * (1/4);
            d_int2_U = secLen * (secNorm.') * derP_U_gphi2(1:2,1:2) * (1/4);
            d_int1_T = secLen * (secNorm.') * derP_T_gphi1(1:2,1:2) * (1/4);
            d_int2_T = secLen * (secNorm.') * derP_T_gphi2(1:2,1:2) * (1/4);
            [ dum1, dum2, dum3, dum4, d_AderB_d_F ] = constLaw( Fe_U, Cp_U, elKU, elGU, etaU, 1, dt, Ame, gphi1e );
            d_int1_U_d_uU = d_AderB_d_F * d_Fe_d_u + d_int1_U * (-1) * (d_uU_p1_d_u + d_uU_p2_d_u); 
            d_int1_U_d_uT = d_int1_U * (d_uT_p1_d_u + d_uT_p2_d_u);
            [ dum1, dum2, dum3, dum4, d_AderB_d_F ] = constLaw( Fe_U, Cp_U, elKU, elGU, etaU, 1, dt, Ame, gphi2e );
            d_int2_U_d_uU = d_AderB_d_F * d_Fe_d_u + d_int2_U * (-1) * (d_uU_p1_d_u + d_uU_p2_d_u); 
            d_int2_U_d_uT = d_int2_U * (d_uT_p1_d_u + d_uT_p2_d_u);
            [ dum1, dum2, dum3, dum4, d_AderB_d_F ] = constLaw( Fe_T, Cp_T, elKT, elGT, etaT, g, dt, Ame, gphi1e );
            d_int1_T_d_uU = d_int1_T * (-1) * (d_uU_p1_d_u + d_uU_p2_d_u);  
            d_int1_T_d_uT = d_AderB_d_F * d_Fe_d_u + d_int1_T * (d_uT_p1_d_u + d_uT_p2_d_u);
            [ dum1, dum2, dum3, dum4, d_AderB_d_F ] = constLaw( Fe_T, Cp_T, elKT, elGT, etaT, g, dt, Ame, gphi2e );
            d_int2_T_d_uU = d_int2_T * (-1) * (d_uU_p1_d_u + d_uU_p2_d_u); 
            d_int2_T_d_uT = d_AderB_d_F * d_Fe_d_u + d_int2_T * (d_uT_p1_d_u + d_uT_p2_d_u);
            TVx_U( (6*jj-5):(6*jj) ) = TVx_U( (6*jj-5):(6*jj) ) - d_int1_U_d_uU;
            TVy_U( (6*jj-5):(6*jj) ) = TVy_U( (6*jj-5):(6*jj) ) - d_int2_U_d_uU;
            TVx_T( (6*jj-5):(6*jj) ) = TVx_T( (6*jj-5):(6*jj) ) - d_int1_T_d_uT;
            TVy_T( (6*jj-5):(6*jj) ) = TVy_T( (6*jj-5):(6*jj) ) - d_int2_T_d_uT;
            TVx_U_ex( (6*jj-5):(6*jj) ) = TVx_U_ex( (6*jj-5):(6*jj) ) - d_int1_U_d_uT;
            TVy_U_ex( (6*jj-5):(6*jj) ) = TVy_U_ex( (6*jj-5):(6*jj) ) - d_int2_U_d_uT;
            TVx_T_ex( (6*jj-5):(6*jj) ) = TVx_T_ex( (6*jj-5):(6*jj) ) - d_int1_T_d_uU;
            TVy_T_ex( (6*jj-5):(6*jj) ) = TVy_T_ex( (6*jj-5):(6*jj) ) - d_int2_T_d_uU;
            
            %. third integral: [phi].[u]
            gam = numPar.NIT_PAR;
            int1 = ( gam/dh ) * secLen * [1 0] * ( 2*phi_p1*uJ_p1 + phi_p1*uJ_p2 + phi_p2*uJ_p1 + 2*phi_p2*uJ_p2 ) * (1/6);
            int2 = ( gam/dh ) * secLen * [0 1] * ( 2*phi_p1*uJ_p1 + phi_p1*uJ_p2 + phi_p2*uJ_p1 + 2*phi_p2*uJ_p2 ) * (1/6);
            %. untransformed
            F_U_od(ii) = F_U_od(ii) - int1;
            F_U_ev(ii) = F_U_ev(ii) - int2;
            %. transformed
            F_T_od(ii) = F_T_od(ii) + int1;
            F_T_ev(ii) = F_T_ev(ii) + int2;
            
            %. third integral, deriv.
            d_int1_d_uU = ( gam/dh ) * secLen * [1 0] * ( 2*phi_p1*(-1)*d_uU_p1_d_u + phi_p1*(-1)*d_uU_p2_d_u + phi_p2*(-1)*d_uU_p1_d_u + 2*phi_p2*(-1)*d_uU_p2_d_u ) * (1/6);
            d_int1_d_uT = ( gam/dh ) * secLen * [1 0] * ( 2*phi_p1*d_uT_p1_d_u + phi_p1*d_uT_p2_d_u + phi_p2*d_uT_p1_d_u + 2*phi_p2*d_uT_p2_d_u ) * (1/6);
            d_int2_d_uU = ( gam/dh ) * secLen * [0 1] * ( 2*phi_p1*(-1)*d_uU_p1_d_u + phi_p1*(-1)*d_uU_p2_d_u + phi_p2*(-1)*d_uU_p1_d_u + 2*phi_p2*(-1)*d_uU_p2_d_u ) * (1/6);
            d_int2_d_uT = ( gam/dh ) * secLen * [0 1] * ( 2*phi_p1*d_uT_p1_d_u + phi_p1*d_uT_p2_d_u + phi_p2*d_uT_p1_d_u + 2*phi_p2*d_uT_p2_d_u ) * (1/6);
            TVx_U( (6*jj-5):(6*jj) ) = TVx_U( (6*jj-5):(6*jj) ) - d_int1_d_uU;
            TVy_U( (6*jj-5):(6*jj) ) = TVy_U( (6*jj-5):(6*jj) ) - d_int2_d_uU;
            TVx_T( (6*jj-5):(6*jj) ) = TVx_T( (6*jj-5):(6*jj) ) + d_int1_d_uT;
            TVy_T( (6*jj-5):(6*jj) ) = TVy_T( (6*jj-5):(6*jj) ) + d_int2_d_uT;
            TVx_U_ex( (6*jj-5):(6*jj) ) = TVx_U_ex( (6*jj-5):(6*jj) ) - d_int1_d_uT;
            TVy_U_ex( (6*jj-5):(6*jj) ) = TVy_U_ex( (6*jj-5):(6*jj) ) - d_int2_d_uT;
            TVx_T_ex( (6*jj-5):(6*jj) ) = TVx_T_ex( (6*jj-5):(6*jj) ) + d_int1_d_uU;
            TVy_T_ex( (6*jj-5):(6*jj) ) = TVy_T_ex( (6*jj-5):(6*jj) ) + d_int2_d_uU;
            
        end
    end

    %. stablisation integrals
    normals = [ 1          0         ;
                0          1         ;
                1/sqrt(2)  1/sqrt(2) ];
    lengths = [ 1 1 2 ];
    for kk=1:size(inde,2)
        
        elInd_loc_R = elem_edeges(inde(kk),1);
        elInd_loc_L = elem_edeges(inde(kk),2);
        elInd_R = elem_inds( elInd_loc_R );
        elInd_L = elem_inds( elInd_loc_L );
        elType_R = elTypes( elInd_R );
        elType_L = elTypes( elInd_L );
        
        if ( elType_R == 3 )||( elType_L == 3 )
            %. one of elements - intersected, stabilise edge
            
            normInd = elem_edeges(inde(kk),3);
            
            omeg = 1;
            if ( elem_edeges(inde(kk),1) < 6.5 )&&( elem_edeges(inde(kk),2) < 6.5 )
                %. node ii belongs to the common edge
                omeg = -1;
            end
            
            nphi1 = (normals(normInd,:).') * [1 0];
            nphi2 = (normals(normInd,:).') * [0 1];
            nphi1e = zeros(3);
            nphi1e(1:2,1:2) = nphi1;
            vnphi1eT = tens2vec(nphi1e.');
            nphi2e = zeros(3);
            nphi2e(1:2,1:2) = nphi2;
            vnphi2eT = tens2vec(nphi2e.');
            
            kap = numPar.STAB_PAR;
            
            d_int1_d_Fe = kap * dh * ( vnphi1eT.' ) * omeg * lengths(normInd);
            d_int2_d_Fe = kap * dh * ( vnphi2eT.' ) * omeg * lengths(normInd);
            
            ind1x = 2*elemNodeInds_all(elInd_R,1)-1;
            ind1y = 2*elemNodeInds_all(elInd_R,1);
            ind2x = 2*elemNodeInds_all(elInd_R,2)-1;
            ind2y = 2*elemNodeInds_all(elInd_R,2);
            ind3x = 2*elemNodeInds_all(elInd_R,3)-1;
            ind3y = 2*elemNodeInds_all(elInd_R,3);
            idx = [ ind1x ind1y ind2x ind2y ind3x ind3y ];
            TR_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) = repmat( ii, 1, 6 );
            TC_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) = idx;
            
            ind1x = 2*elemNodeInds_all(elInd_L,1)-1;
            ind1y = 2*elemNodeInds_all(elInd_L,1);
            ind2x = 2*elemNodeInds_all(elInd_L,2)-1;
            ind2y = 2*elemNodeInds_all(elInd_L,2);
            ind3x = 2*elemNodeInds_all(elInd_L,3)-1;
            ind3y = 2*elemNodeInds_all(elInd_L,3);
            idx = [ ind1x ind1y ind2x ind2y ind3x ind3y ];
            TR_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) = repmat( ii, 1, 6 );
            TC_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) = idx;

            elemNodeInds_R = elemNodeInds_all( elInd_R, : );
            d_Fe_d_u_R = elemNodeInds_R(4) * d_Fe_d_u_r * (1/dh);
            elemNodeInds_L = elemNodeInds_all( elInd_L, : );
            d_Fe_d_u_L = elemNodeInds_L(4) * d_Fe_d_u_r * (1/dh);
            
            if ( elType_R ~= 1 )&&( elType_L ~= 1 )
                %. both elements not transformed, stabilise untransformed
                
                Fe_U_R = vec2tens( Fe_U_all( elInd_R, : ) );
                Fe_U_L = vec2tens( Fe_U_all( elInd_L, : ) );
                
                jumpF_U = Fe_U_R - Fe_U_L;
                
                int1_U = kap * dh * tensTrace2( nphi1e * jumpF_U ) * omeg * lengths(normInd);
                int2_U = kap * dh * tensTrace2( nphi2e * jumpF_U ) * omeg * lengths(normInd);
                
                F_U_od(ii) = F_U_od(ii) + int1_U;
                F_U_ev(ii) = F_U_ev(ii) + int2_U;
                
                TVx_U_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) = TVx_U_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) + d_int1_d_Fe * d_Fe_d_u_R;
                TVy_U_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) = TVy_U_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) + d_int2_d_Fe * d_Fe_d_u_R;
                
                TVx_U_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) = TVx_U_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) - d_int1_d_Fe * d_Fe_d_u_L;
                TVy_U_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) = TVy_U_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) - d_int2_d_Fe * d_Fe_d_u_L;
                
            end
            if ( elType_R ~= 2 )&&( elType_L ~= 2 )
                %. both elements not untransformed, stabilise transformed
                
                Fe_T_R = vec2tens( Fe_T_all( elInd_R, : ) );
                Fe_T_L = vec2tens( Fe_T_all( elInd_L, : ) );
                
                jumpF_T = Fe_T_R - Fe_T_L;
                
                int1_T = kap * dh * tensTrace2( nphi1e * jumpF_T ) * omeg * lengths(normInd);
                int2_T = kap * dh * tensTrace2( nphi2e * jumpF_T ) * omeg * lengths(normInd);
                
                F_T_od(ii) = F_T_od(ii) + int1_T;
                F_T_ev(ii) = F_T_ev(ii) + int2_T;
                
                TVx_T_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) = TVx_T_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) + d_int1_d_Fe * d_Fe_d_u_R;
                TVy_T_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) = TVy_T_n( (6*elInd_loc_R-5):(6*elInd_loc_R) ) + d_int2_d_Fe * d_Fe_d_u_R;
                
                TVx_T_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) = TVx_T_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) - d_int1_d_Fe * d_Fe_d_u_L;
                TVy_T_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) = TVy_T_n( (6*elInd_loc_L-5):(6*elInd_loc_L) ) - d_int2_d_Fe * d_Fe_d_u_L;
                
            end
        end
    end
    
    JR(ii, :) = TR;
    JC(ii, :) = TC;
    JVx_U(ii, :) = TVx_U;
    JVy_U(ii, :) = TVy_U;
    JVx_T(ii, :) = TVx_T;
    JVy_T(ii, :) = TVy_T;
    JVx_U_ex(ii, :) = TVx_U_ex;
    JVy_U_ex(ii, :) = TVy_U_ex;
    JVx_T_ex(ii, :) = TVx_T_ex;
    JVy_T_ex(ii, :) = TVy_T_ex;
    JR_n(ii, :) = TR_n;
    JC_n(ii, :) = TC_n;
    JVx_U_n(ii, :) = TVx_U_n;
    JVy_U_n(ii, :) = TVy_U_n;
    JVx_T_n(ii, :) = TVx_T_n;
    JVy_T_n(ii, :) = TVy_T_n;
    
end

F_U_all( 2:2:(2*Nx*Ny), : ) = F_U_ev;
F_U_all( 1:2:(2*Nx*Ny-1), : ) = F_U_od;
F_T_all( 2:2:(2*Nx*Ny), : ) = F_T_ev;
F_T_all( 1:2:(2*Nx*Ny-1), : ) = F_T_od;
F = [ F_U_all; F_T_all; ];

JRr = JR(:);
JCr = JC(:);
JVx_Ur = JVx_U(:);
JVy_Ur = JVy_U(:);
JVx_Tr = JVx_T(:);
JVy_Tr = JVy_T(:);
JVx_U_exr = JVx_U_ex(:);
JVy_U_exr = JVy_U_ex(:);
JVx_T_exr = JVx_T_ex(:);
JVy_T_exr = JVy_T_ex(:);
idxrem = find( JRr==0 );
JRr( idxrem ) = [];
JCr( idxrem ) = [];
JVx_Ur( idxrem ) = [];
JVy_Ur( idxrem ) = [];
JVx_Tr( idxrem ) = [];
JVy_Tr( idxrem ) = [];
JVx_U_exr( idxrem ) = [];
JVy_U_exr( idxrem ) = [];
JVx_T_exr( idxrem ) = [];
JVy_T_exr( idxrem ) = [];
JR_nr = JR_n(:);
JC_nr = JC_n(:);
JVx_U_nr = JVx_U_n(:);
JVy_U_nr = JVy_U_n(:);
JVx_T_nr = JVx_T_n(:);
JVy_T_nr = JVy_T_n(:);
idxrem = find( JR_nr==0 );
JR_nr( idxrem ) = [];
JC_nr( idxrem ) = [];
JVx_U_nr( idxrem ) = [];
JVy_U_nr( idxrem ) = [];
JVx_T_nr( idxrem ) = [];
JVy_T_nr( idxrem ) = [];
RR = [ 2*JRr-1; 2*JRr; 2*Nx*Ny+2*JRr-1; 2*Nx*Ny+2*JRr; 2*JRr-1; 2*JRr; 2*Nx*Ny+2*JRr-1; 2*Nx*Ny+2*JRr; 2*JR_nr-1; 2*JR_nr; 2*Nx*Ny+2*JR_nr-1; 2*Nx*Ny+2*JR_nr; ];
CC = [ JCr; JCr; 2*Nx*Ny+JCr; 2*Nx*Ny+JCr; 2*Nx*Ny+JCr; 2*Nx*Ny+JCr; JCr; JCr; JC_nr; JC_nr; 2*Nx*Ny+JC_nr; 2*Nx*Ny+JC_nr; ];
VV = [ JVx_Ur; JVy_Ur; JVx_Tr; JVy_Tr; JVx_U_exr; JVy_U_exr; JVx_T_exr; JVy_T_exr; JVx_U_nr; JVy_U_nr; JVx_T_nr; JVy_T_nr; ];
Js = sparse( RR, CC, VV );

%% prescribe unused nodes

F(2*exclNodes_U-1) = u(2*exclNodes_U-1);
F(2*exclNodes_U) = u(2*exclNodes_U);

F(2*Nx*Ny+2*exclNodes_T-1) = u(2*Nx*Ny+2*exclNodes_T-1);
F(2*Nx*Ny+2*exclNodes_T) = u(2*Nx*Ny+2*exclNodes_T);

Js( 2*exclNodes_U-1, : ) = 0;
Js( sub2ind( size(Js), 2*exclNodes_U-1, 2*exclNodes_U-1 ) ) = 1;
Js( 2*exclNodes_U, : ) = 0;
Js( sub2ind( size(Js), 2*exclNodes_U, 2*exclNodes_U ) ) = 1;

Js( 2*Nx*Ny+2*exclNodes_T-1, : ) = 0;
Js( sub2ind( size(Js), 2*Nx*Ny+2*exclNodes_T-1, 2*Nx*Ny+2*exclNodes_T-1 ) ) = 1;
Js( 2*Nx*Ny+2*exclNodes_T, : ) = 0;
Js( sub2ind( size(Js), 2*Nx*Ny+2*exclNodes_T, 2*Nx*Ny+2*exclNodes_T ) ) = 1;

%% prescribe coinciding nodes

F(2*intNodesUnk_C-1) = F(2*intNodesUnk_C-1) + F(2*Nx*Ny+2*intNodesUnk_C-1);
F(2*intNodesUnk_C) = F(2*intNodesUnk_C) + F(2*Nx*Ny+2*intNodesUnk_C);
Js( 2*intNodesUnk_C-1, : ) = Js( 2*intNodesUnk_C-1, : ) + Js( 2*Nx*Ny+2*intNodesUnk_C-1, : );
Js( 2*intNodesUnk_C, : ) = Js( 2*intNodesUnk_C, : ) + Js( 2*Nx*Ny+2*intNodesUnk_C, : );

F(2*Nx*Ny+2*intNodesUnk_C-1) = u(2*Nx*Ny+2*intNodesUnk_C-1) - u(2*intNodesUnk_C-1);
F(2*Nx*Ny+2*intNodesUnk_C) = u(2*Nx*Ny+2*intNodesUnk_C) - u(2*intNodesUnk_C);

Js( 2*Nx*Ny+2*intNodesUnk_C-1, : ) = 0;
Js( sub2ind( size(Js), 2*Nx*Ny+2*intNodesUnk_C-1, 2*Nx*Ny+2*intNodesUnk_C-1 ) ) = 1;
Js( sub2ind( size(Js), 2*Nx*Ny+2*intNodesUnk_C-1, 2*intNodesUnk_C-1 ) ) = -1;
Js( 2*Nx*Ny+2*intNodesUnk_C, : ) = 0;
Js( sub2ind( size(Js), 2*Nx*Ny+2*intNodesUnk_C, 2*Nx*Ny+2*intNodesUnk_C ) ) = 1;
Js( sub2ind( size(Js), 2*Nx*Ny+2*intNodesUnk_C, 2*intNodesUnk_C ) ) = -1;

%% boudary conditions

indLx_U = (2:(2*Nx):(2*Nx*(Ny-1)+2))-1;
indLy_U = (2:(2*Nx):(2*Nx*(Ny-1)+2));
indRx_U = (2:(2*Nx):(2*Nx*(Ny-1)+2))+2*Nx-3;
indRy_U = (2:(2*Nx):(2*Nx*(Ny-1)+2))+2*Nx-2;
indBx_U = (2:2:(2*Nx))-1;
indBy_U = (2:2:(2*Nx));
indTx_U = (2:2:(2*Nx))+2*Nx*(Ny-1)-1;
indTy_U = (2:2:(2*Nx))+2*Nx*(Ny-1);

indLBx_U = 1;
indLBy_U = 2;
indRBx_U = 2*Nx-1;
indRBy_U = 2*Nx;
indLTx_U = 2*Nx*(Ny-1)+1;
indLTy_U = 2*Nx*(Ny-1)+2;
indRTx_U = 2*Nx*Ny-1;
indRTy_U = 2*Nx*Ny;

indLx_T = (2:(2*Nx):(2*Nx*(Ny-1)+2))+2*Nx*Ny-1;
indLy_T = (2:(2*Nx):(2*Nx*(Ny-1)+2))+2*Nx*Ny;
indRx_T = (2:(2*Nx):(2*Nx*(Ny-1)+2))+2*Nx*Ny+2*Nx-3;
indRy_T = (2:(2*Nx):(2*Nx*(Ny-1)+2))+2*Nx*Ny+2*Nx-2;
indBx_T = (2:2:(2*Nx))+2*Nx*Ny-1;
indBy_T = (2:2:(2*Nx))+2*Nx*Ny;
indTx_T = (2:2:(2*Nx))+2*Nx*Ny+2*Nx*(Ny-1)-1;
indTy_T = (2:2:(2*Nx))+2*Nx*Ny+2*Nx*(Ny-1);

indLBx_T = 2*Nx*Ny+1;
indLBy_T = 2*Nx*Ny+2;
indRBx_T = 2*Nx*Ny+2*Nx-1;
indRBy_T = 2*Nx*Ny+2*Nx;
indLTx_T = 2*Nx*Ny+2*Nx*(Ny-1)+1;
indLTy_T = 2*Nx*Ny+2*Nx*(Ny-1)+2;
indRTx_T = 2*Nx*Ny+2*Nx*Ny-1;
indRTy_T = 2*Nx*Ny+2*Nx*Ny;

xx = dh*(0:1:(Nx-1)).';
yy = dh*(0:1:(Ny-1)).';

if ( geomPar.BCtype == 0 )
    %. no rigid body motion, enforced at BL and BR corners (transf.)

    F(indLBx_T) = u(indLBx_T);
    F(indLBy_T) = u(indLBy_T);
    F(indRBy_T) = u(indRBy_T);

    Js( indLBx_T, : ) = 0;
    Js( indLBx_T, indLBx_T ) = 1;
    Js( indLBy_T, : ) = 0;
    Js( indLBy_T, indLBy_T ) = 1;
    Js( indRBy_T, : ) = 0;
    Js( indRBy_T, indRBy_T ) = 1;

elseif ( geomPar.BCtype == 1 )
    %. clamped bottom (transf.) and inclined indentation top (untransf.)

    %. B.C. transf. bottom
    F(indBx_T) = u(indBx_T);
    F(indBy_T) = u(indBy_T);

    %. B.C. untransf. top
    F(indTx_U) = u(indTx_U);
    F(indTy_U) = u(indTy_U) + 0.1*xx;

    Js( indBx_T, : ) = 0;
    Js( sub2ind( size(Js), indBx_T, indBx_T ) ) = 1;
    Js( indBy_T, : ) = 0;
    Js( sub2ind( size(Js), indBy_T, indBy_T ) ) = 1;
    Js( indTx_U, : ) = 0;
    Js( sub2ind( size(Js), indTx_U, indTx_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;

elseif ( geomPar.BCtype == 2 )
    %. clamped bottom (untransf.) and inclined indentation top (untransf.)

    %. B.C. untransf. bottom
    F(indBx_U) = u(indBx_U);
    F(indBy_U) = u(indBy_U);

    %. B.C. untransf. top
    F(indTx_U) = u(indTx_U);
    F(indTy_U) = u(indTy_U) + 0.1*xx;

    Js( indBx_U, : ) = 0;
    Js( sub2ind( size(Js), indBx_U, indBx_U ) ) = 1;
    Js( indBy_U, : ) = 0;
    Js( sub2ind( size(Js), indBy_U, indBy_U ) ) = 1;
    Js( indTx_U, : ) = 0;
    Js( sub2ind( size(Js), indTx_U, indTx_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;

elseif ( geomPar.BCtype == 3 )
    %. sliding bottom/left/right/top, biaxial stretching (transf.)

    displ = geomPar.BCpar_displ;
    
    %. B.C. transf. bottom
    F(indBy_T) = u(indBy_T);

    %. B.C. transf. left
    F(indLx_T) = u(indLx_T);
    
    %. B.C. transf. right
    F(indRx_T) = u(indRx_T) - displ;
    
    %. B.C. transf. top
    F(indTy_T) = u(indTy_T) - displ;

    Js( indBy_T, : ) = 0;
    Js( sub2ind( size(Js), indBy_T, indBy_T ) ) = 1;
    Js( indLx_T, : ) = 0;
    Js( sub2ind( size(Js), indLx_T, indLx_T ) ) = 1;
    Js( indRx_T, : ) = 0;
    Js( sub2ind( size(Js), indRx_T, indRx_T ) ) = 1;
    Js( indTy_T, : ) = 0;
    Js( sub2ind( size(Js), indTy_T, indTy_T ) ) = 1;

elseif ( geomPar.BCtype == 4 )
    %. no rigid body motion, enforced at BL and BR corners (untransf.)

    F(indLBx_U) = u(indLBx_U);
    F(indLBy_U) = u(indLBy_U);
    F(indRBy_U) = u(indRBy_U);

    Js( indLBx_U, : ) = 0;
    Js( indLBx_U, indLBx_U ) = 1;
    Js( indLBy_U, : ) = 0;
    Js( indLBy_U, indLBy_U ) = 1;
    Js( indRBy_U, : ) = 0;
    Js( indRBy_U, indRBy_U ) = 1;
    
elseif ( geomPar.BCtype == 5 )
    %. complex boundary displacement (transf.)

    a = geomPar.BCpar_a;
    b = geomPar.BCpar_b;
    c = geomPar.BCpar_c;
    d = geomPar.BCpar_d;
    
    L2 = dh*(Nx-1)/2;
    H2 = dh*(Ny-1)/2;
    uxT = (xx-L2)*b/L2;
    uyT = a + abs(xx-L2)*(c-a)/L2;
    uxR = d - abs(yy-H2)*(d-b)/H2;
    uyR = (yy-H2)*c/H2;
    
    %. B.C. transf. bottom
    F(indBx_T) = u(indBx_T) - uxT;
    F(indBy_T) = u(indBy_T) + uyT;

    %. B.C. transf. left
    F(indLx_T) = u(indLx_T) + uxR;
    F(indLy_T) = u(indLy_T) - uyR;
    
    %. B.C. transf. right
    F(indRx_T) = u(indRx_T) - uxR;
    F(indRy_T) = u(indRy_T) - uyR;
    
    %. B.C. transf. top
    F(indTx_T) = u(indTx_T) - uxT;
    F(indTy_T) = u(indTy_T) - uyT;

    Js( indBx_T, : ) = 0;
    Js( sub2ind( size(Js), indBx_T, indBx_T ) ) = 1;
    Js( indBy_T, : ) = 0;
    Js( sub2ind( size(Js), indBy_T, indBy_T ) ) = 1;
    Js( indLx_T, : ) = 0;
    Js( sub2ind( size(Js), indLx_T, indLx_T ) ) = 1;
    Js( indLy_T, : ) = 0;
    Js( sub2ind( size(Js), indLy_T, indLy_T ) ) = 1;
    Js( indRx_T, : ) = 0;
    Js( sub2ind( size(Js), indRx_T, indRx_T ) ) = 1;
    Js( indRy_T, : ) = 0;
    Js( sub2ind( size(Js), indRy_T, indRy_T ) ) = 1;
    Js( indTx_T, : ) = 0;
    Js( sub2ind( size(Js), indTx_T, indTx_T ) ) = 1;
    Js( indTy_T, : ) = 0;
    Js( sub2ind( size(Js), indTy_T, indTy_T ) ) = 1;

elseif ( geomPar.BCtype == 6 )
    %. clamped bottom (transf.), sliding left/right (transf.), sliding left/right (untransf.), clamped and displaced top (untransf.)

    displ = geomPar.BCpar_displ;
    
    %. B.C. transf. bottom
    F(indBx_T) = u(indBx_T);
    F(indBy_T) = u(indBy_T);

    %. B.C. transf. left
    F(indLx_T) = u(indLx_T);
    
    %. B.C. transf. right
    F(indRx_T) = u(indRx_T);
    
    %. B.C. untransf. left
    F(indLx_U) = u(indLx_U);
    
    %. B.C. untransf. right
    F(indRx_U) = u(indRx_U);
    
    %. B.C. untransf. top
    F(indTx_U) = u(indTx_U);
    F(indTy_U) = u(indTy_U) - displ;

    Js( indBx_T, : ) = 0;
    Js( sub2ind( size(Js), indBx_T, indBx_T ) ) = 1;
    Js( indBy_T, : ) = 0;
    Js( sub2ind( size(Js), indBy_T, indBy_T ) ) = 1;
    Js( indLx_T, : ) = 0;
    Js( sub2ind( size(Js), indLx_T, indLx_T ) ) = 1;
    Js( indRx_T, : ) = 0;
    Js( sub2ind( size(Js), indRx_T, indRx_T ) ) = 1;
    Js( indLx_U, : ) = 0;
    Js( sub2ind( size(Js), indLx_U, indLx_U ) ) = 1;
    Js( indRx_U, : ) = 0;
    Js( sub2ind( size(Js), indRx_U, indRx_U ) ) = 1;
    Js( indTx_U, : ) = 0;
    Js( sub2ind( size(Js), indTx_U, indTx_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;

elseif ( geomPar.BCtype == 9 )
    %. clamped bottom (transf.) and displaced top (untransf.)

    displ_h = geomPar.BCpar_displ_h;
    displ_v = geomPar.BCpar_displ_v;
    
    %. B.C. transf. bottom
    F(indBx_T) = u(indBx_T);
    F(indBy_T) = u(indBy_T);

    %. B.C. untransf. top
    F(indTx_U) = u(indTx_U) - displ_h;
    F(indTy_U) = u(indTy_U) - displ_v;

    Js( indBx_T, : ) = 0;
    Js( sub2ind( size(Js), indBx_T, indBx_T ) ) = 1;
    Js( indBy_T, : ) = 0;
    Js( sub2ind( size(Js), indBy_T, indBy_T ) ) = 1;
    Js( indTx_U, : ) = 0;
    Js( sub2ind( size(Js), indTx_U, indTx_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;

elseif ( geomPar.BCtype == 10 )
    %. sliding bottom (untransf.) and displaced sliding top (untransf.)

    displ = geomPar.BCpar_displ;

    %. B.C. untransf. BL corner
    F(indLBx_U) = u(indLBx_U);

    %. B.C. untransf. bottom
    F(indBy_U) = u(indBy_U);

    %. B.C. untransf. top
    F(indTy_U) = u(indTy_U) - displ;

    Js( indLBx_U, : ) = 0;
    Js( indLBx_U, indLBx_U ) = 1;
    Js( indBy_U, : ) = 0;
    Js( sub2ind( size(Js), indBy_U, indBy_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;

elseif ( geomPar.BCtype == 11 )
    %. sliding bottom/left/right/top, constrained "uniaxial" stretching (untransf.)

    displ = geomPar.BCpar_displ;
    
    %. B.C. transf. bottom
    F(indBy_U) = u(indBy_U);

    %. B.C. transf. left
    F(indLx_U) = u(indLx_U);
    
    %. B.C. transf. right
    F(indRx_U) = u(indRx_U) + displ/(1+displ);
    
    %. B.C. transf. top
    F(indTy_U) = u(indTy_U) - displ;

    Js( indBy_U, : ) = 0;
    Js( sub2ind( size(Js), indBy_U, indBy_U ) ) = 1;
    Js( indLx_U, : ) = 0;
    Js( sub2ind( size(Js), indLx_U, indLx_U ) ) = 1;
    Js( indRx_U, : ) = 0;
    Js( sub2ind( size(Js), indRx_U, indRx_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;

elseif ( geomPar.BCtype == 15 )
    %. patch test B.C., untransf., for unit square geom.

    %. B
    F(indBx_U) = u(indBx_U) + xx*0.05;
    F(indBy_U) = u(indBy_U) - xx*0.1;
    
    %. L
    F(indLx_U) = u(indLx_U) - yy*0.02;
    F(indLy_U) = u(indLy_U) + yy*0.01;
    
    %. R
    F(indRx_U) = u(indRx_U) - yy*0.02 + 0.05;
    F(indRy_U) = u(indRy_U) + yy*0.01 - 0.1;
    
    %. T
    F(indTx_U) = u(indTx_U) + xx*0.05 - 0.02;
    F(indTy_U) = u(indTy_U) - xx*0.1 + 0.01;
    
    Js( indBx_U, : ) = 0;
    Js( sub2ind( size(Js), indBx_U, indBx_U ) ) = 1;
    Js( indBy_U, : ) = 0;
    Js( sub2ind( size(Js), indBy_U, indBy_U ) ) = 1;
    Js( indLx_U, : ) = 0;
    Js( sub2ind( size(Js), indLx_U, indLx_U ) ) = 1;
    Js( indLy_U, : ) = 0;
    Js( sub2ind( size(Js), indLy_U, indLy_U ) ) = 1;
    Js( indRx_U, : ) = 0;
    Js( sub2ind( size(Js), indRx_U, indRx_U ) ) = 1;
    Js( indRy_U, : ) = 0;
    Js( sub2ind( size(Js), indRy_U, indRy_U ) ) = 1;
    Js( indTx_U, : ) = 0;
    Js( sub2ind( size(Js), indTx_U, indTx_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;
    
elseif ( geomPar.BCtype == 16 )
    %. clamped bottom (transf.), sliding left/right (transf.), sliding left/right (untransf.), clamped and displaced top with a ramp (untransf.)

    displ = geomPar.BCpar_displ;
    rm = geomPar.BCpar_rm;

    L = dh*(Nx-1);
    uy = displ + xx*rm/L;
    
    %. B.C. transf. bottom
    F(indBx_T) = u(indBx_T);
    F(indBy_T) = u(indBy_T);

    %. B.C. transf. left
    F(indLx_T) = u(indLx_T);
    
    %. B.C. transf. right
    F(indRx_T) = u(indRx_T);
    
    %. B.C. untransf. left
    F(indLx_U) = u(indLx_U);
    
    %. B.C. untransf. right
    F(indRx_U) = u(indRx_U);
    
    %. B.C. untransf. top
    F(indTx_U) = u(indTx_U);
    F(indTy_U) = u(indTy_U) - uy;

    Js( indBx_T, : ) = 0;
    Js( sub2ind( size(Js), indBx_T, indBx_T ) ) = 1;
    Js( indBy_T, : ) = 0;
    Js( sub2ind( size(Js), indBy_T, indBy_T ) ) = 1;
    Js( indLx_T, : ) = 0;
    Js( sub2ind( size(Js), indLx_T, indLx_T ) ) = 1;
    Js( indRx_T, : ) = 0;
    Js( sub2ind( size(Js), indRx_T, indRx_T ) ) = 1;
    Js( indLx_U, : ) = 0;
    Js( sub2ind( size(Js), indLx_U, indLx_U ) ) = 1;
    Js( indRx_U, : ) = 0;
    Js( sub2ind( size(Js), indRx_U, indRx_U ) ) = 1;
    Js( indTx_U, : ) = 0;
    Js( sub2ind( size(Js), indTx_U, indTx_U ) ) = 1;
    Js( indTy_U, : ) = 0;
    Js( sub2ind( size(Js), indTy_U, indTy_U ) ) = 1;

end

end

