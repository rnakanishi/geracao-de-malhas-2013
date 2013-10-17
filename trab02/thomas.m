%
% Em parceria com Lucas Pagliosa
%
function [Err, h]= thomas()

    [m, n, Rb, Rt, Rl, Rr] = domain();

	% Discretizando valores xi e eta
	xi = linspace(0.,1,m) ;
	eta = linspace(0.,1.,n) ;
	DXi = xi(2)-xi(1);
	DEta = eta(2)-eta(1);

	% Alocando matrizes para os eixos X e Y
	X = zeros(m,n) ;
	Y = zeros(m,n) ;

	% TRANSFINITE INTERPOLATION
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

	[a,b,c,d,lineXi, lineEta, pointXi, pointEta] = ttmparams();

	P = zeros(m, n);
	Q = zeros(m, n);
	for i = 2:m-1
	  for j = 2:n-1
	    P(i, j) = control(xi(i), eta(j), lineXi, pointXi, pointEta, a, b, c, d);
	    Q(i, j) = control(eta(j), xi(i), lineEta, pointEta, pointXi, a, b, c, d);
	  end
	end
	display(P);
	display(Q);
	
	%to declare
	err = inf;
	tol = 1e-4;
	iter = 0;	

	p =[]; qx=[]; qy=[];
	Err = [];
    
	while (err > tol & iter < 1000)


		iter = iter+1;
		err = 0;
		for j=2:n-1
			for i=2:m-1
				dxXi = X(i+1,j)-X(i-1,j);
				dyXi = Y(i+1,j)-Y(i-1,j);
				dxEta = X(i,j+1)-X(i,j-1);
				dyEta = Y(i,j+1)-Y(i,j-1);

				g11 = (dxXi^2 + dyXi^2)/(4*DXi^2);
				g22 = (dxEta^2 + dyEta^2)/(4*DEta^2);
				g12 = (1/4*DXi*DEta)*(dxXi*dxEta + dyXi*dyEta);
				g = g11*g22 - g12^2;

				d2xXi = dxXi^2/(DXi^2);
				d2xEta = dxEta^2/DEta^2;
				d2xXiEta = (X(i+1,j+1)+X(i-1,j-1)-X(i-1,j+1)-X(i+1,j-1))/(4*DXi*DEta);
				d2yXi = dyXi^2/(DXi^2);
				d2yEta = dyEta^2/DEta^2;
				d2yXiEta = (Y(i+1,j+1)+Y(i-1,j-1)-Y(i-1,j+1)-Y(i+1,j-1))/(4*DXi*DEta);

				a = g22/(DXi^2) - g*P(i,j)/(2*DXi);
				c = g22/(DXi^2) + g*P(i,j)/(2*DXi);
				b = 2*g22/(DXi^2) + 2*g11/(DEta^2);
				d = g11*(X(i+1,j)+X(i-1,j))/(DEta^2) - 2*g12*d2xXiEta - g*Q(i,j)*(dxEta);
				e = g11*(Y(i,j+1)+Y(i,j-1))/(DEta^2) - 2*g12*d2yXiEta - g*Q(i,j)*(dyEta);


				if i==2
					p(i) = c/b ;
					qx(i) = (a*X(i-1,j) + d)/b;
					qy(i) = (a*Y(i-1,j) + e)/b;
				elseif i==m-1
					p(i) = 0;
					qx(i) = (c*X(i+1,j) + d + a*qx(i-1))/(b - a*p(i-1));
					qy(i) = (c*Y(i+1,j) + e + a*qy(i-1))/(b - a*p(i-1));
				else
					p(i) = c/(b - a*p(i-1));
					qx(i) = ( d + a*qx(i-1))/(b - a*p(i-1));
					qy(i) = ( e + a*qy(i-1))/(b - a*p(i-1));
				end
			end

			for it = m-1:-1:2
				xtemp(it,j) = p(it)*X(it+1,j) + qx(it);
				ytemp(it,j) = p(it)*Y(it+1,j) + qy(it);

				err = err + hypot( (X(it,j) - xtemp(it,j)), (Y(it,j) - ytemp(it,j)) );

				X(it,j) = xtemp(it,j);
				Y(it,j) = ytemp(it,j);
			end			% err = sqrt(err);
		end	
		Err = [Err err];
	end
	fprintf(1,'Numero de iteracoes de Thomas: %d\n',iter);
	% Plotando
	h=figure;
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

function xnew = solveThomas(a,b,c,d,x)
	n = length(d);

	c(1) = c(1) / (b(1)+1); %Divisao por zero
	d(1) = d(1) / (b(1)+1);
	for i=2:n-1
		div = b(i) - a(i-1)*c(i-1);
		c(i) = c(i) / div;
		d(i) = (d(i) - a(i-1)*d(i-1))/div;
	end
	d(n) = (d(n) - a(n-1)*d(n-1))/(b(n) - a(n-1)*c(n-1));

	x(n+1) = d(n);
	for i=n-1:-1:1
		x(i+1) = d(i) - c(i) * x(i+2);
	end
	xnew = x;
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
        % fprintf(1,'point it\n');
        value = value - sign( xi - pointXi(point) )*b(point)*...
            exp(-d(point)*hypot( xi - pointXi(point),eta - pointEta(point)));
    end
    % fprintf(1,'%d\n', value);
end

