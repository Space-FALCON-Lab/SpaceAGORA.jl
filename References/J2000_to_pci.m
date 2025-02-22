% VSCode readable version of the code in the J2000_to_pci.mlx file in the
% References folder. This script is used to derive the transformation matrix
% from the J2000 ECI frame to the planet-centric inertial frame for an
% arbitrary planet given the right ascension and declination of the planet's
% North pole of rotation. The notebook file will be more readable in MATLAB, 
% but this script is provided for those who prefer to read the code in a text editor.

syms ra dec
earth_north = [0;0;1]; % North pole of rotation of the Earth in the J2000 ECI frame

% North pole of rotation of an arbitrary planet given right ascension and declination
planet_north = [cos(dec)*cos(ra); cos(dec)*sin(ra); sin(dec)]; 

% The first axis in the planet-centric inertial frame is the axis along the
% intersection of the planet's equatorial plane with the Earth's equatorial
% plane, which can be defined as the cross product of the Earth's North
% pole of rotation and the planet's North pole of rotation
i_hat = cross(earth_north, planet_north);
i_hat = i_hat/norm(i_hat); % Normalize the result of the cross product

% The third axis is just the planet's North pole of rotation
k_hat = planet_north; 

% The second axis is derived from the cross product of the third axis,
% which is just the planet's North pole of rotation vector, and the first
% axis defined above
j_hat = cross(k_hat, i_hat);
j_hat = j_hat/norm(j_hat);

pci_rot = simplify([i_hat j_hat k_hat])' % Transpose to get orientation matrix from rotation matrix

% The expression derived here can be simplified further, which has been done in the 
% J2000_to_pci matrix defined in Planet_data.jl, but I'm not good enough with MATLAB 
% to get it to show the fully simplified version. Essentially, sigma_2 was
% eliminated since it just simplifies to cos(dec), which then cancels with
% the same term in the numerator everywhere it's used, and sigma_1 was
% simplified by combining the cos(ra) and sin(ra) terms.

