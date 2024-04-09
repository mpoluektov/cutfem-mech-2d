function [ We_U_all_av, We_T_all_av, Fe_U_all_av, Fe_T_all_av, Pe_U_all_av, Pe_T_all_av ] = interElemAver( resElem, geomPar, intParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nx = geomPar.Nx;
Ny = geomPar.Ny;

We_U_all = resElem.We_U_all;
We_T_all = resElem.We_T_all;
Fe_U_all = resElem.Fe_U_all;
Fe_T_all = resElem.Fe_T_all;
Pe_U_all = resElem.Pe_U_all;
Pe_T_all = resElem.Pe_T_all;

We_U_all_av = zeros( Nx*Ny, 1 );
We_T_all_av = zeros( Nx*Ny, 1 );
Fe_U_all_av = zeros( Nx*Ny, 9 );
Fe_T_all_av = zeros( Nx*Ny, 9 );
Pe_U_all_av = zeros( Nx*Ny, 9 );
Pe_T_all_av = zeros( Nx*Ny, 9 );

intElem = intParam.intElem;
fracElem = intParam.fracElem;
elTypes = intParam.elTypes;

for ii=1:(Nx*Ny)
               
    %. squares
    sq_ar = ii - floor( (ii-0.5)/Nx );
    sq_al = ii - floor( (ii-0.5)/Nx ) - 1;
    sq_br = ii - floor( (ii-0.5)/Nx ) - (Nx-1);
    sq_bl = ii - floor( (ii-0.5)/Nx ) - (Nx-1) - 1;

    elem_inds = [ 2*sq_al-1 ;
                  2*sq_ar-1 ;
                  2*sq_br-1 ;
                    2*sq_al ;
                    2*sq_br ;
                    2*sq_bl ];

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

    area_U = zeros( size(inda,2), 1 );
    area_T = zeros( size(inda,2), 1 );
    We_U = zeros( size(inda,2), 1 );
    We_T = zeros( size(inda,2), 1 );
    Fe_U = zeros( size(inda,2), 9 );
    Fe_T = zeros( size(inda,2), 9 );
    Pe_U = zeros( size(inda,2), 9 );
    Pe_T = zeros( size(inda,2), 9 );
    
    for jj=1:size(inda,2)
        
        elType_loc = elTypes( elem_inds(inda(jj)) );
        if ( elType_loc == 2 )
            %. untransformed
            area_U(jj) = 1;
            area_T(jj) = 0;
        elseif ( elType_loc == 1 )
            %. transformed
            area_U(jj) = 0;
            area_T(jj) = 1;
        else
            %. intersected
            useFracArea = 1;
            if ( useFracArea == 1 )
                indEl_loc = find( intElem == elem_inds(inda(jj)), 1, 'first' );
                area_U(jj) = fracElem(indEl_loc);
                area_T(jj) = 1 - fracElem(indEl_loc);
            else
                area_U(jj) = 1;
                area_T(jj) = 1;
            end
        end
        
        We_U(jj,:) = We_U_all( elem_inds(inda(jj)), : );
        We_T(jj,:) = We_T_all( elem_inds(inda(jj)), : );
        Fe_U(jj,:) = Fe_U_all( elem_inds(inda(jj)), : );
        Fe_T(jj,:) = Fe_T_all( elem_inds(inda(jj)), : );
        Pe_U(jj,:) = Pe_U_all( elem_inds(inda(jj)), : );
        Pe_T(jj,:) = Pe_T_all( elem_inds(inda(jj)), : );
        
    end
    
    area_U_sum = sum(area_U);
    area_T_sum = sum(area_T);
    if ( area_U_sum ~= 0 )
        We_U_all_av(ii,:) = area_U.' * We_U / area_U_sum;
        Fe_U_all_av(ii,:) = area_U.' * Fe_U / area_U_sum;
        Pe_U_all_av(ii,:) = area_U.' * Pe_U / area_U_sum;
    end
    if ( area_T_sum ~= 0 )
        We_T_all_av(ii,:) = area_T.' * We_T / area_T_sum;
        Fe_T_all_av(ii,:) = area_T.' * Fe_T / area_T_sum;
        Pe_T_all_av(ii,:) = area_T.' * Pe_T / area_T_sum;
    end
end

end

