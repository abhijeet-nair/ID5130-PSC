function [u,v] = getIndVel(g,x,y,xo,yo)
    % To get induced velocities due to a given vortex, placed at (xo,yo) 
    % at a point (x,y)
    % Input: 
    % g - Vortex strength
    % x,y - Location at which velocity is estimated
    % xo,yo - Location at which vortex is present
    % Output:
    % u - x direction velocity
    % v - y direction velocity
   
    den = (x-xo).^2 + (y-yo).^2;
    u = (g/(2*pi)).*(y-yo)./den;
    v = -(g/(2*pi)).*(x-xo)./den;
end
