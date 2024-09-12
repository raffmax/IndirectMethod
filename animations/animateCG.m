function frames = animateCG(X, footout, gamma, filename)
% animateCG Compute the hybrid flow of a given model.
%   FRAMES = ANIMATE(T, X0) returns animated FRAMES starting from X0
%   until time t based on the flow map, jump map, jump set, and drawing
%   functions.
%
%   [FRAMES, SOL] = ANIMATE(___) returns FRAMES and the ode45 output
%   structure SOL.
%
% See also: ODE45, FLOW


%% Animate
% We draw frames to the current axes and save them to a video file

% animation parameters
nframes = size(X,1); % # of frames to draw between [0, t] units of time

ax = gca; % autoscales by default
axis equal
grid off
set(ax,'XTickLabel',[]);
set(ax,'YTickLabel',[]);


% save frames throughout the animation
frames(nframes) = struct('cdata',[],'colormap',[]); 

% here we get a graphics object with no parent
graphic = draw();
set(graphic, 'Parent', ax);

% draw initial state
draw(X(1,:), gamma, [0; 0], graphic);
frames(1) = getframe(ax);

% animation loop
for i = 2:nframes
    % update drawing
    foot = footout(:, i);
    if ~isempty(foot)
        foot = foot(:, end);
    else
        foot = [0; 0];
    end
    
    draw(X(i,:)', gamma, foot, graphic);
    drawnow limitrate;
    % save drawing
    frames(i) = getframe(ax);
end


if ~isempty(filename)
    movie2gif(frames, [filename,'.gif'],'DelayTime',.05,'LoopCount',Inf)
end
end