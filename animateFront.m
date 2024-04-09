
clear variables;
close all;

%. palette
col1 = [0,160,176]/255;
col2 = [106,74,60]/255;
col3 = [204,51,63]/255;
col4 = [235,104,65]/255;
col5 = [237,201,81]/255;
col6 = [0.7,0.7,0.7];
col7 = [0,0,0];

mycolmap = [0,0,0;
            165,0,38;
            215,48,39;
            244,109,67;
            253,174,97;
            254,224,144;
            255,255,191;
            224,243,248;
            171,217,233;
            116,173,209;
            69,117,180;
            49,54,149;
            127,127,127]/255;
mycolmap = mycolmap(end:(-1):1,:);

%. symbols
symb_eps = char(949);
symb_thet = char(952);
symb_pi = char(960);
symb_sig = char(963);
symb_phi = char(966);
symb_psi = char(968);
symb_omega = char(969);
symb_Gamma = char(915);

fps = 25;
total_t = 4;

%. video output
fn = 'tmp1';
mov = VideoWriter( fn, 'MPEG-4' );
mov.FrameRate = fps;
mov.Quality = 100;
open(mov);

%. data input
fn = 'tmp1';
load( [ fn '/param' ] );

t_min = 0;
t_max = tend;
Nframes = round(total_t*fps);

scale = 2.5;

aspRatio = 1;
fontname = 'Segoe UI';

axisXlen = 60; %. mm
axisXlen_pt = axisXlen / 0.352778; %. pt
tickLen = 1; %. mm
tickLen_norm = tickLen / axisXlen; %. frac.
    
figure( 'Units', 'points', 'Position', [ 100 40 160*scale+2*axisXlen_pt*scale 70*scale+axisXlen_pt*scale/aspRatio ] );

set( gcf, 'Color', 'w' );

Nx = geomPar.Nx;
Ny = geomPar.Ny;

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

t_all = 0:dt:tend;

scaleX_min = 0;
scaleX_max = 1;
scaleY_min = 0;
scaleY_max = 1;

for ii = 0:Nframes
    
    clf;

    subplot(1,2,1);
    ax1 = gca;
    subplot(1,2,2);
    ax2 = gca;
    
    t_cur = t_min + (t_max-t_min) * ii/Nframes;
    [ dum, ind ] = min( abs( t_all - t_cur ) );
    load( [ fn '/res' sprintf('%05i',ind) ], 'res' );
    
    %%
    
    axes(ax1);
    hold on;

    triplot( tri, coordGrid(1,:)', coordGrid(2,:)', '-', 'LineWidth', 0.5*scale, 'Color', col6 );
    
    intP = res.intParam.intPoints;
    for jj=1:size(res.intParam.intElem,1)
        n1 = res.intParam.intConn(jj,1);
        n2 = res.intParam.intConn(jj,2);
        plot( [ intP(n1,1) intP(n2,1) ], [ intP(n1,2) intP(n2,2) ], '-', 'LineWidth', 0.5*scale, 'Color', col3 );
    end

    text( 0, scaleY_max+0.05, 'reference configuration', 'FontName', fontname, 'FontSize', 10*scale, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom' );
    text( 0, scaleY_max+0.15, [ '{\itt} = ' sprintf( '%.1f', t_all(ind) ) ], 'FontName', fontname, 'FontSize', 10*scale, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom' );
    
    hold off;
    
    pbaspect( [ aspRatio 1 1 ] );
    limits = [ scaleX_min scaleX_max scaleY_min scaleY_max ];
    axis( limits );
    axis off;

    set( gca, 'Units', 'points' );
    set( gca, 'Position', [ 50*scale 20*scale axisXlen_pt*scale axisXlen_pt*scale/aspRatio ] );
    
    %%
    
    axes(ax2);
    hold on;

    elTypes = res.intParam.elTypes;
    We = zeros(2*(Nx-1)*(Ny-1),1);
    We(elTypes==1,:) = res.We_T_elem(elTypes==1);
    We(elTypes==2,:) = res.We_U_elem(elTypes==2);
    We(elTypes==3,:) = ( res.We_T_elem(elTypes==3) + res.We_U_elem(elTypes==3) )/2;
    
    vals = We;
    limMin = 0;
    limMax = 0.06;
    Ncm = size(mycolmap,1)-2;
    dcol = (limMax-limMin)/Ncm;
    valsSc = round( 1+(Ncm-1)*(vals-limMin-dcol/2)/(limMax-limMin-dcol) );
    valsSc = valsSc+1;
    valsSc( valsSc>(Ncm+1) ) = Ncm+2;
    valsSc( valsSc<2 ) = 1;
    trisurf( tri, coordGrid(1,:)', coordGrid(2,:)', zeros(Nx*Ny,1), valsSc, 'CDataMapping', 'direct', 'LineStyle', 'none' );
    colormap(ax2,mycolmap);
    pos = [ 110*scale+2*axisXlen_pt*scale 20*scale 5*scale axisXlen_pt*scale/aspRatio ];
    Nticks = Ncm+3;
    ticks = [ 2, (Nticks+1)/2, Nticks-1 ];
    labs = { sprintf('%.2f',limMin) sprintf('%.2f',(limMin+limMax)/2) sprintf('%.2f',limMax) };
    colorbar( 'LineWidth', 0.5*scale, 'FontName', fontname, 'FontSize', 10*scale, 'Box', 'off', 'TickDirection', 'out', 'Units', 'points', 'Position', pos, 'Ticks', ticks, 'TickLabels', labs );

    text( 0, scaleY_max+0.05, 'strain energy density', 'FontName', fontname, 'FontSize', 10*scale, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom' );
    
    hold off;
    
    pbaspect( [ aspRatio 1 1 ] );
    limits = [ scaleX_min scaleX_max scaleY_min scaleY_max ];
    axis( limits );
    axis off;

    set( gca, 'Units', 'points' );
    set( gca, 'Position', [ 100*scale+axisXlen_pt*scale 20*scale axisXlen_pt*scale axisXlen_pt*scale/aspRatio ] );
    
    %%
    
    F = getframe(gcf);
    writeVideo( mov, F );
end

close(mov);

