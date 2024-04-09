function [ ] = mainPhaseTran( filename )
%MAINPHASETRAN Cut finite elements for modelling phase transitions
%   See (Poluektov and Figiel, Comput. Mech., 63:885-911, 2019).
%   Code has been written by M. Poluektov.

%% param.

if ( nargin < 1 )
    error( 'Specify filename' );
end

%. elastic properties
matPar.elGU = 24;
matPar.elGT = 10;
matPar.elKU = 75;
matPar.elKT = 32;

%. viscosities, not used in elasticity
matPar.etaU = 1;
matPar.etaT = 1;

%. chemical parameters
matPar.alpD = 2;
matPar.ksD = 0.1;
matPar.gam = -0.15;
matPar.coefANN = 0.018;
matPar.coefVel = 0.043;

%. phase transition (0) or chemical reaction (1)
matPar.isChem = 0;

%. grid
domainSize = 1;
geomPar.Nx = 24 + 1;
geomPar.Ny = 24 + 1;
geomPar.dh = domainSize/(geomPar.Nx-1);

%. boundary conditions for mechanics
geomPar.BCtype = 3;
geomPar.BCpar_displ = 0.038;

%. boundary conditions for diffusion
geomPar.BCdiff = 1;

%. numerical parameters
numPar.ATOL = 1e-11;
numPar.DTOL = 1e-2;
numPar.STOL = 1e-2;
numPar.MAX_ITER = 30;
numPar.SDLT = 1e-7;
numPar.NIT_PAR = 1e4;
numPar.STAB_PAR = 1e-1;
numPar.PTOL = 1e-14;

%. time
tend = 100;
Nt = 100 + 1;
dt = tend/(Nt-1);

%. transformation ratio
g_loadcases_ref = 1.05;

%% initialise

Nx = geomPar.Nx;
Ny = geomPar.Ny;
dh = geomPar.dh;

%. create grid
gXs = ( 0:1:(Nx-1) )*dh;
gYs = ( 0:1:(Ny-1) )*dh;
gX = repmat(gXs,1,Ny);
gYa = repmat(gYs,1,Nx);
gYb = reshape(gYa,Ny,[]);
gYt = gYb.';
gY = gYt(:).';
coordGrid = [ gX; gY; ];

%% create initial interface

intType = 1;

if ( intType == 0 )
    %. flat interface
    
    geomPar.windThresh = 0;
    
    y0 = 0.1;
    
    intP = [  -dh y0 ;
             1+dh y0 ];

    conn = [ 1 2 ];

elseif ( intType == 1 )
    %. circular interface

    geomPar.windThresh = 0.5;
    
    x0 = 0.5;
    y0 = 0.5;
    a = 0.27;

    dalp = pi/180;
    st = 0:360;
    px = x0 + a*cos(dalp*st);
    py = y0 + a*sin(dalp*st);
    intP = [ px.' py.' ];
    
    NintP = size( intP, 1 );
    conn = [ 1:(NintP-1) ;
             2:NintP     ].';

end

ppStart = 1;

%% plot

plotInt = 0;
if ( plotInt == 1 )
    tri = zeros( 2*(Nx-1)*(Ny-1), 3 );
    for ii=1:((Nx-1)*(Ny-1))

        ind_bl = ii + floor((ii-0.5)/(Nx-1));
        ind_br = ii + 1 + floor((ii-0.5)/(Nx-1));
        ind_tl = ii + Nx + floor((ii-0.5)/(Nx-1));
        ind_tr = ii + Nx + 1 + floor((ii-0.5)/(Nx-1));

        tri( 2*ii-1, 1 ) = ind_br;
        tri( 2*ii-1, 2 ) = ind_bl;
        tri( 2*ii-1, 3 ) = ind_tl;
        tri( 2*ii, 1 ) = ind_tl;
        tri( 2*ii, 2 ) = ind_tr;
        tri( 2*ii, 3 ) = ind_br;
    end

    figure(1);
    hold on;

    triplot( tri, coordGrid(1,:)', coordGrid(2,:)', 'Color', [0.75 0.75 0.75] );

    plot( intP(:,1), intP(:,2), 'd-' );

    limits = [ -0.1 1.1 -0.1 1.1 ];
    axis( limits );
    pbaspect( [ 1 1 1 ] );
end

%% init. plasticity

Cp_prev_U_all = zeros( 2*(Nx-1)*(Ny-1), 9 );
Cp_prev_T_all = zeros( 2*(Nx-1)*(Ny-1), 9 );
Cp_prev_U_all(:,1:3) = 1;
Cp_prev_T_all(:,1:3) = 1;
prevState.Cp_prev_U_all = Cp_prev_U_all;
prevState.Cp_prev_T_all = Cp_prev_T_all;

%% solve and move front

mkdir( filename );
save( [ filename '/param' ], 'coordGrid', 'matPar', 'geomPar', 'numPar', 'Nt', 'tend', 'dt', 'g_loadcases_ref' );

for pp=ppStart:Nt

    %. intersection points and elements
    intParam = calcIntLevelSet( intP, conn, geomPar, numPar.PTOL );

    if ( plotInt == 1 )
        %. transf. elem.
        triplot( tri(intParam.elTypes==1,:), coordGrid(1,:)', coordGrid(2,:)', 'Color', 'b' );
        %. untransf. elem.
        triplot( tri(intParam.elTypes==2,:), coordGrid(1,:)', coordGrid(2,:)', 'Color', 'm' );
    end
    
    %. initial estimate
    if ( pp == 1 )
        u0 = zeros( 4*Nx*Ny, 1 );
        g_loadcases = g_loadcases_ref;
    else
        %. better initial estimate
        urUt = ucmr;
        urTt = ucmr;
        urUt(:,intParam.exclNodes_U) = 0;
        urTt(:,intParam.exclNodes_T) = 0;
        ur = [ urUt urTt ];
        u0 = ur(:);
        g_loadcases = g_loadcases_ref(end);
    end

    %. solve mechanical
    [ u, noConv_mech, resElem ] = solveMechPhaseTran( u0, matPar, geomPar, numPar, intParam, prevState, g_loadcases, dt );
    if ( noConv_mech ~= 0 )
        break;
    end
    
    %. reshape solution
    ur = reshape(u,2,[]);
    urUt = ur(:,1:(Nx*Ny));
    urTt = ur(:,(Nx*Ny+1):(2*Nx*Ny));
    urUt(:,intParam.exclNodes_U) = 0;
    urUt(:,intParam.intNodesUnk_U) = 0;
    urUt(:,intParam.intNodesUnk_C) = urUt(:,intParam.intNodesUnk_C)/2;
    urTt(:,intParam.exclNodes_T) = 0;
    urTt(:,intParam.intNodesUnk_T) = 0;
    urTt(:,intParam.intNodesUnk_C) = urTt(:,intParam.intNodesUnk_C)/2;
    ucmr = urUt + urTt;
    
    %. displacement at the interface
    NintP = size(intParam.intPoints,1);
    u_int = zeros(NintP,2);
    for ii=1:NintP
        n1 = intParam.intNodes(ii,1);
        n2 = intParam.intNodes(ii,2);
        w1 = intParam.fracFaces(ii);
        w2 = 1 - intParam.fracFaces(ii);
        un1 = ucmr(:,n1).';
        un2 = ucmr(:,n2).';
        u_int(ii,:) = un1*w1 + un2*w2 ;
    end
    
    %. energy at the interface
    [ We_U_all_av, We_T_all_av, Fe_U_all_av, Fe_T_all_av, Pe_U_all_av, Pe_T_all_av ] = interElemAver( resElem, geomPar, intParam );
    [ WUi, WTi, PjumpFi ] = intEnerg( We_U_all_av, We_T_all_av, Fe_U_all_av, Fe_T_all_av, Pe_U_all_av, Pe_T_all_av, intParam );
    xi = matPar.gam - WTi + WUi + PjumpFi;

    %. solve diffusion
    if ( matPar.isChem == 1 )
        %. chemo-mechanics
        ceq = exp( -matPar.coefANN*xi );
        conc = solveDiff( ceq, matPar, geomPar, numPar, intParam );
    else
        %. phase transition
        ceq = 0;
        conc = 0;
    end

    %. move interface
    intP_new = moveInt( conc, ceq, xi, matPar, intParam, dt );

    if ( plotInt == 1 )
        plot( intP_new(:,1), intP_new(:,2), 'd' );
    end

    %. export
    res.u = u;
    res.ucmr = ucmr;
    res.u_int = u_int;
    res.We_U_elem = resElem.We_U_all;
    res.We_T_elem = resElem.We_T_all;
    res.Fe_U_elem = resElem.Fe_U_all;
    res.Fe_T_elem = resElem.Fe_T_all;
    res.Pe_U_elem = resElem.Pe_U_all;
    res.Pe_T_elem = resElem.Pe_T_all;
    res.Cp_U_elem = resElem.Cp_U_all;
    res.Cp_T_elem = resElem.Cp_T_all;
    res.We_U_node = We_U_all_av;
    res.We_T_node = We_T_all_av;
    res.Fe_U_node = Fe_U_all_av;
    res.Fe_T_node = Fe_T_all_av;
    res.Pe_U_node = Pe_U_all_av;
    res.Pe_T_node = Pe_T_all_av;
    res.intParam = intParam;
    res.WUi = WUi;
    res.WTi = WTi;
    res.PjumpFi = PjumpFi;
    res.conc = conc;
    res.intP_new = intP_new;

    %. update interface points
    intP = intP_new;
    conn = intParam.intConn;
    
    %. update plastic deformations
    prevState.Cp_prev_U_all = resElem.Cp_U_all;
    prevState.Cp_prev_T_all = resElem.Cp_T_all;
    elTypes = intParam.elTypes;
    indU = (elTypes==2);
    indT = (elTypes==1);
    prevState.Cp_prev_U_all(indT,:) = resElem.Cp_T_all(indT,:);
    prevState.Cp_prev_T_all(indU,:) = resElem.Cp_U_all(indU,:);
    
    fprintf( '  step %.0f/%.0f done\n', pp, Nt );
    
    %. save
    save( [ filename '/res' sprintf('%05i',pp) ], 'res' );

end

end

