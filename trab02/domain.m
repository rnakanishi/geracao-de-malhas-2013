function [m, n, Rb, Rt, Rl, Rr] = domain
    % Resolucao dos eixos xi e eta
    m = 20 ;
    n = 20 ;

	% Fronteiras do dominio fisico

	% Waves
	% Rb = @(s) ([s, 15-17*s-1.3*s^2]);
	% Rt = @(s) ([s, 13+14*s+2*s^2]);
	% Rl = @(s) ([-1-.1*sin(s*7*pi/2), s]);
	% Rr = @(s) ([1+.1*sin(s*7*pi/2), s]);

	% M
	Rb = @(s) ([s,-0.2*cos(s*4*pi)]);
	Rt = @(s) ([s,0.2*cos(s* 2*pi)]);
	Rl = @(s) ([s-1.7, s]);
	Rr = @(s) ([-s+1.7, s]);

	% FISH
	% Rb = @(s) ([s,-abs(0.2-s)*double(0|s<0.2) + double(0|s>=0.2)*(-sin((s-0.2)*(pi-pi/8)))]);
	% Rt = @(s) ([s,abs(0.2-s)*double(0|s<0.2) + double(0|s>=0.2)*(sin((s-0.2)*(pi-pi/8)))]);,
	% Rl = @(s) ([0, s]);
	% Rr = @(s) ([abs(s-0.5)+1, s]);
end