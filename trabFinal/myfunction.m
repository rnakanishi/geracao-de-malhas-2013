function [U,S] = myfunction(filename)

    if nargin<1
        display('Usage: myfunction( <mesh_file> );');
        return;
    end

    options.name = filename;
    [vertex,faces] = read_mesh(filename);
    n = size(vertex,2);
    
    options.symmetrize = 1;
    options.normalize = 0;
    display('Computing Laplacian Matrix');
    L = compute_mesh_laplacian(vertex,faces,'conformal',options);

    if size(L,1) > 1e4
        display ( ['Laplacian matrix too big: ' int2str(size(L,1))] );
        return
    end
    display('Computing eigenvectors');
    [U,S] = eig(full(L)); S = diag(S);
    display('sorting');
    [S,I] = sort(S,'ascend'); U = U(:,I);
    
    fiedler = U(:,2);
    
    options.face_vertex_color = fiedler;

    h = figure;
    set(h,'name', filename(8:length(filename)-4),'numbertitle','off');
    colormap jet;
    set(h, 'Position',[0 0 1500 500]);
    clf;
    subplot(1,3,1);
    plot_mesh(vertex,faces,options);
        
    % Property of mesh
    [~, I] = sort(fiedler,'ascend');
    [umin,umax,cmin,cmax,cmean,cgauss,normal] = compute_curvature(vertex,faces);
    property = cmean(I);
    
%     figure;
    % Scattergrams
    NN = size(L,1);
    lf = ceil(log2(NN));
    
    subplot(1,3,2);
    plot(1:NN,property);
    title('Curvatura Média');
    pbaspect([1,1,1]);
    
    y = zeros(1,2^(lf));
    y(1:NN) = property;
    
    N = 2^(lf);
    y = y';
       
    display('Calculating wavelet filter');
    window_size = 32;
    filt_opt = default_filter_options('audio', window_size);
 
    % Only compute zeroth-, first- and second-order scattering.
    order = 2;
    scat_opt.M = order;

    % Prepare wavelet transforms to use in scattering.
    [Wop, ~] = wavelet_factory_1d(N, filt_opt, scat_opt);
%     figure;
    % Compute the scattering coefficients of y.
    display('Applying scattering');
    S = scat(y, Wop);

    % Display first-order coefficients, as a function of time and scale j1, and
    % second-order coefficients for a fixed j1, as a function of time and second
    % scale j2.
    j1 = 23;
    img = scattergram(S{2},[]);
    img = img{1};    
    
    imagesc(img(:,1:ceil(NN/(window_size/2))+1));
    title('Scattergram');
    pbaspect([1,1,1]);
%     % Renormalize second order by dividing it by the first order and compute the
%     % logarithm of the coefficients.
%     S = renorm_scat(S);
%     S = log_scat(S);
% 
%     % Display the transformed coefficients.
%     scattergram(S{2},[],S{3},j1);
    
    colormap jet;
%     saveas(h,['../figures/' filename(8:length(filename)-3) 'png'],'png');
end