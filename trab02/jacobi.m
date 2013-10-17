function [Err, h]=  jacobi

    [m, n, Rb, Rt, Rl, Rr] = domain();


    % Discretizando valores xi e eta
    xi = linspace(0.,1,m) ;
    eta = linspace(0.,1.,n) ;
    

    % Alocando matrizes para os eixos X e Y
    X = zeros(m,n) ;
    Y = zeros(m,n) ;

    % Fronteiras do dominio fisico

    
    for i = 1:m
        Xi = xi(i) ;
        for j = 1:n
            Eta = eta(j) ;
            
            % TFI
            XY = (1-Eta)*Rb(Xi)+Eta*Rt(Xi)...
               + (1-Xi)*Rl(Eta)+Xi*Rr(Eta)...;
               - (Xi*Eta*Rt(1)+Xi*(1-Eta)*Rb(1)+Eta*(1-Xi)*Rt(0)+(1-Xi)*(1-Eta)*Rb(0)) ;
        
            X(i,j) = XY(1) ;
            Y(i,j) = XY(2) ;
            
        end
    end


    % Equacoes de Winslow
    k=0;
    a=zeros(m-2,n-2);b=a;c=a;
    err = inf;
    tol = 1e-4;
    iter = 0;
    % figure;
    [a,b,c,d,lineXi, lineEta, pointXi, pointEta] = ttmparams();
    
    P = zeros(m, n);
    Q = zeros(m, n);
    for i = 2:m-1
      for j = 2:n-1
        P(i, j) = control(xi(i), eta(j), lineXi, pointXi, pointEta, a, b, c, d);
        Q(i, j) = control(eta(j), xi(i), lineEta, pointEta, pointXi, a, b, c, d);
      end
    end
    % display(P);
    % display(Q);

    Err = [];
    xtemp = X; ytemp = Y;
    while (err > tol)
        err = 0;
        iter = iter + 1;
        for j=2:n-1
            for i=2:m-1
                
                pp = P(i,j);
                qq = Q(i,j);

                % Tensores metricos
                c(i-1,j-1)=((X(i+1,j)-X(i-1,j))^2 + (Y(i+1,j)-Y(i-1,j))^2) ; %g11
                a(i-1,j-1)=((X(i,j+1)-X(i,j-1))^2 + (Y(i,j+1)-Y(i,j-1))^2) ; %g22
                b(i-1,j-1)=(X(i+1,j)-X(i-1,j)) * (X(i,j+1)-X(i,j-1)) + (Y(i+1,j)-Y(i-1,j)) * (Y(i,j+1)-Y(i,j-1)); %g12
                g(i-1,j-1)= a(i-1,j-1)*c(i-1,j-1) - b(i-1,j-1)^2;

                xtemp(i,j) = 1/(2*(a(i-1,j-1)+c(i-1,j-1)))*(...
                    a(i-1,j-1)*X(i+1,j) - 0.5*b(i-1,j-1)*X(i+1,j+1) + 0.5*b(i-1,j-1)*X(i+1,j-1)+ ...
                    c(i-1,j-1)*X(i,j+1) + c(i-1,j-1)*X(i,j-1) + ...
                    a(i-1,j-1)*X(i-1,j) - 0.5*b(i-1,j-1)*X(i-1,j-1) + 0.5*b(i-1,j-1)*X(i-1,j+1) + ...
                    g(i-1,j-1)*( pp*(X(i+1,j) -X(i-1,j)) + (qq*(X(i,j+1)-X(i,j-1)) ) ) )...
                    ;
                
                ytemp(i,j) = 1/(2*(a(i-1,j-1)+c(i-1,j-1)))*(...
                    a(i-1,j-1)*Y(i+1,j) - 0.5*b(i-1,j-1)*Y(i+1,j+1) + 0.5*b(i-1,j-1)*Y(i+1,j-1)+ ...
                    c(i-1,j-1)*Y(i,j+1) + c(i-1,j-1)*Y(i,j-1) + ...
                    a(i-1,j-1)*Y(i-1,j) - 0.5*b(i-1,j-1)*Y(i-1,j-1) + 0.5*b(i-1,j-1)*Y(i-1,j+1) + ...
                    g(i-1,j-1)*( pp*(Y(i+1,j)-Y(i-1,j)) + (qq*(Y(i,j+1)+Y(i,j-1)) ) ) )...
                    ;
                              
          end
       end
        err = (norm(X-xtemp)+norm(Y-ytemp));
        X = xtemp;
        Y = ytemp;
        Err = [Err err];
       %  clc;display(err);
       % pause(0.001);
       % clf
       % hold on
       %  set(gcf,'color','w') ;
       %  axis equal;
       %  axis off
       %  box on
       %  for i=1:m
       %      plot(X(i,:),Y(i,:),'k','linewidth',1);
       %  end
       %  for j=1:n
       %      plot(X(:,j),Y(:,j),'k','linewidth',1);
       %  end
       %  hold off
       %  M(iter) = getframe;
       % err = sqrt(err);
    end

    fprintf(1,'Numero de iteracoes de Jacobi: %d\n',iter);
    % Plotando
    h = figure;
    clf
    set(gcf,'color','w') ;
    axis equal;
    axis off
    box on
    hold on
    for i=1:m
        plot(X(i,:),Y(i,:),'k','linewidth',1);
    end
    for j=1:n
        plot(X(:,j),Y(:,j),'k','linewidth',1);
    end
    hold off
   
end

function value = control(xi, eta, lineXi, pointXi, pointEta, a, b, c, d)
    
    value = 0;
    % display(lineXi);
    % XI
    % fprintf(1,'line it, %d lines\n', length(lineXi));
    for line = 1:length(lineXi)
        value = value - ((sign( xi - lineXi(line) ).*a(line)).*...
            exp(-c(line) .*abs( xi - lineXi(line) )));
    end
    size_ = size(pointXi);
    for point = 1:size_(1)
        value = value - sign( xi - pointXi(point) )*b(point)*...
            exp(-d(point)*hypot( xi - pointXi(point),eta - pointEta(point)));
    end
    % fprintf(1,'%d\n', value);
end

