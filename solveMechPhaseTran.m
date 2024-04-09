function [ u, noConv, resElem ] = solveMechPhaseTran( u0, matPar, geomPar, numPar, intParam, prevState, g_loadcases, dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Nx = geomPar.Nx;
Ny = geomPar.Ny;

ATOL = numPar.ATOL;
DTOL = numPar.DTOL;
STOL = numPar.STOL;
MAX_ITER = numPar.MAX_ITER;
SDLT = numPar.SDLT;

g_ind = 1;
g_cur = g_loadcases(g_ind);
g_lcSize = size(g_loadcases,2);

u = u0;

noConv = 1;
nrIter = 1;
while ( noConv == 1 )
    
    [ F, J_an, resElem ] = calcMechResPhaseTran( u, matPar, geomPar, numPar, intParam, prevState, g_cur, dt );
    
    %. numerical Jacobian
    calcNumJ = 0;
    if ( calcNumJ == 1 )
        J_num = zeros(4*Nx*Ny,4*Nx*Ny);
        for kk = 1:(4*Nx*Ny)
            u_c = u;
            u_c(kk) = u_c(kk) + SDLT;
            [ F_c ] = calcMechResPhaseTran( u_c, matPar, geomPar, numPar, intParam, prevState, g_cur, dt );
            J_num( :, kk ) = ( F_c - F ) / SDLT;
        end
        J_num_s = sparse(J_num);
    end

    F_s = sparse(F);

    Dlt_s = -J_an \ F_s;
    Dlt = full(Dlt_s);
    
    normF = max(abs(F));
    normDlt = max(abs(Dlt));
    
    regSol = 0;
    lamb = 1;
    if ( max(normF,normDlt) > DTOL )&&( regSol == 1 )
        lamb = 0.5;
    end

    u = u + lamb*Dlt;
    
    if ( normF < ATOL )&&( normDlt < ATOL )&&( g_ind == g_lcSize )&&( isreal(Dlt) )
        noConv = 0;
        fprintf( '    N.R. converged, %.0f iter., normF=%1.1e, normDlt=%1.1e\n', nrIter, normF, normDlt );
    else
        fprintf( '        iteration completed, normF=%1.1e, normDlt=%1.1e, lamb=%0.1f, isreal(Dlt)=%0.0f, loadcase=%0.0f/%0.0f \n', normF, normDlt, lamb, isreal(Dlt), g_ind, g_lcSize );
    end
    if ( nrIter > MAX_ITER )
        noConv = 2;
        fprintf( '    NO CONVERGENCE, N.R.\n' );
    end
    if ( ~isreal(Dlt) )
        noConv = 3;
        fprintf( '    NO CONVERGENCE, N.R.\n' );
    end
    nrIter = nrIter + 1;
    
    if ( max(normF,normDlt) < STOL )&&( g_ind < g_lcSize )
        g_ind = g_ind + 1;
        g_cur = g_loadcases(g_ind);
    end
    
end

end

