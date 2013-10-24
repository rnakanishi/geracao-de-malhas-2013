function [U,S] = myfunction(filename)
    
    options.name = filename;
    [vertex,faces] = read_off(filename);
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
    
    options.face_vertex_color = U(:,2);

    h = figure;
    set(h, 'Position',[30 30 1500 720]);
    clf;
    subplot(2,3,[1,4]);
    plot_mesh(vertex,faces,options);
    
    N = size(L,1);
    lf = floor(log2(N));
    
%     figure;
    subplot(2,3,[2,5]);
    plot(1:N,U(:,2));
    
    y = zeros(1,2^(lf+1));
    y(1:N) = U(:,2);
    
    N = 2^(lf+1);
    y = y';
        
    filt_opt = default_filter_options('audio', 32);
 
    % Only compute zeroth-, first- and second-order scattering.
    order = 2;
    scat_opt.M = order;

    % Prepare wavelet transforms to use in scattering.
    [Wop, ~] = wavelet_factory_1d(N, filt_opt, scat_opt);
%     figure;
    % Compute the scattering coefficients of y.
    S = scat(y, Wop);

    % Display first-order coefficients, as a function of time and scale j1, and
    % second-order coefficients for a fixed j1, as a function of time and second
    % scale j2.
    j1 = 23;
    scattergram(S{2},[],S{3},j1);

    % Renormalize second order by dividing it by the first order and compute the
    % logarithm of the coefficients.
    S = renorm_scat(S);
    S = log_scat(S);

    % Display the transformed coefficients.
    scattergram(S{2},[],S{3},j1);
end