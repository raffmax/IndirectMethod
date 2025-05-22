function h = draw(x, gamma, foot, h)
% DRAW Draws the rimless wheel.
% H = DRAW() draws and returns a graphical handle of the rimless wheel at
% the wheel's home position (i.e., x = 0) and a graphical handle of the
% rolling surface.
%
% H = DRAW(X, H) updates and returns the graphic handle H of the rimless
% wheel and rolling surface at state X.  If there exists an axes object
% that is an ancestor of H, then the rolling surface is updated to fit the
% axes limits.
%
% The rimless wheel consists of a hub and spokes that rolls on a surface.

if nargin == 0
    h = new_walker();
else    
    p = parameters();
    L = p.L;
    % current angles
    theta = x(1);
    phi = x(2);  
        
    % convert angles into drawing coordinate system
    a = [-gamma; theta; -phi]; 
    
    hip = [-L*sin(theta - gamma); +L*sin(theta - gamma)*tan(gamma)];
    if foot(1) >0
        %hip = [-L*sin(-theta+phi - gamma); +L*sin(theta-phi - gamma)*tan(gamma)];
        hip = hip+[0; -L*cos(phi - theta)+L*cos(theta)];
    end
        
    % get bipd, stance, and swing
    biped = findobj(h(2), 'Tag', 'biped');
    stance = findobj(h(2), 'Tag', 'stance');
    swing = findobj(h(2), 'Tag', 'swing');
    
    % rotate about origin and tranlsate to actual pivot point
    T = makehgtform('translate', [-hip' 0], 'zrotate', a(1));
    set(biped, 'Matrix', T);
    
    T = makehgtform('zrotate', a(2));
    set(stance, 'Matrix', T);
        
    T = makehgtform('translate', [0 L 0], 'zrotate', a(3));
    set(swing, 'Matrix', T);
    
    % update surface
    draw_ground(a(1), h(1));
end
end

function h = new_walker()
% NEW_WALKER Draws a new graphic of the rimless wheel at its home position.
% H = NEW_WALKER() returns a graphics array consisting of a graphical object
% representing the wheel as the first element and the rolling surface as
% the second element in the array H.
%
% See also DRAW_GROUND

%% Model Parameters and Colors
p = parameters();
L = p.L; % scaled leg length
r_hip = L / 8; % hip radius
r_foot = r_hip / 3; % foot radius
w = r_foot;

% https://www.beschaeftigte.uni-stuttgart.de/uni-services/
%    oeffentlichkeitsarbeit/corporate-design/cd-dateien/manual_englisch.pdf
stance_color = [0,0.2549,0.5686];%'#004191'; % midblue
swing_color =  [0.1961,0.1961,0.1961];%'#323232'; % anthracite
mass_color =    [0,0.7451,1];% '#00BEFF'; % lightblue

%% The Biped
% We draw the biped as an hgtransform that translates and rotates the
% subgraphics (hip and legs) attached to it.  This includes rotating the
% wheel so it is on the correct slope.

%% biped
biped = hgtransform('Parent', [], 'Tag', 'biped');

%% stance leg
% we'll draw legs at origin so that rotations are correct with respect to
% model definition and then translate into correct position
stance = hgtransform('Parent', biped, 'Tag', 'stance');
pos = [-w/4 0 w/2 L]';
stance_leg = rectangle('Parent', stance, 'Position', pos, 'Curvature', 0.2);
set(stance_leg, 'FaceColor', stance_color);

pos = [-1 -1 2 2] * r_foot + [0 L/2 0 0];
stance_foot = rectangle('Parent', stance, 'Position', pos, 'Curvature', [1 1]);
set(stance_foot, 'FaceColor', mass_color);

pos = [-1 -1 2 2] * r_hip + [0 L 0 0];
hip = rectangle('Parent', stance, 'Position', pos, 'Curvature', [1 1]);
set(hip, 'FaceColor', mass_color);

%% swing leg
T = makehgtform('translate', [0 L 0]);
swing = hgtransform('Parent', stance, 'Matrix', T, 'Tag', 'swing');

pos = [-w/4 -L w/2 L]';
swing_leg = rectangle('Parent', swing, 'Position', pos, 'Curvature', 0.2);
set(swing_leg, 'FaceColor', swing_color);

pos = [-1 -1 2 2] * r_foot + [0 -L/2 0 0];
swing_foot = rectangle('Parent', swing, 'Position', pos, 'Curvature', [1 1]);
set(swing_foot, 'FaceColor', mass_color);

%% surface
% The surface is a line, which we initially draw as a line segment that is
% the length of the wheel.  On future updates, we draw the line to extend
% to the entire axes limits using draw_ground.

ground = line('Parent', [], 'LineWidth', 1.5, ...
    'XData', [-L, L], 'YData', [0, 0]);

h = [ground,biped];
end