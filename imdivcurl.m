function [divergence, curl] = imdivcurl(x, y, u, v)
% Computes the divergence and curl of a 2D velocity field.
% x, y: 2D arrays of x and y coordinates
% u, v: 2D arrays of x and y velocity components
% divergence: 2D array of velocity field divergence
% curl: 2D array of velocity field curl

% Compute the x and y derivatives of the velocity components
u_x = gradient(u, x(1,:), y(:,1));
v_x = gradient(v, x(1,:), y(:,1));
u_y = gradient(u, y(:,1), x(1,:));
v_y = gradient(v, y(:,1), x(1,:));

% Compute the velocity field divergence and curl
divergence = u_x + v_y;
curl = v_x - u_y;

end

% this a test