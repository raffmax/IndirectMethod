function h = draw_ground(gamma, h)
% GROUND Draws the ground at incline gamma.
% H = GROUND() returns a line object with zero width.
%
% H = GROUND(GAMMA, H) returns a line object with incline gamma visible in
% the limits of an axes object (set as an ancestor of H).

if nargin == 0
    h = line('Parent', [], 'LineWidth', 1.5, ...
        'XData', [-0, 0], 'YData', [0, 0]);
else
    % update surface
    ax = ancestor(h, 'axes');
    if ~isempty(ax)
        % draw "infinite line"
        [X, Y] = grounddata(gamma, ax);
        set(h, 'XData', X, 'YData', Y);
    end
end
end

function [X, Y] = grounddata(gamma, ax)
% GROUNDDATA returns limits for line object to fit inside axes limit.
% [X, Y] = GROUNDDATE(AX, GAMMA) returns the end points X and Y of a line
% passing through the origin with slope GAMMA in the visible window defined
% by the axes limits of AX.

% get line info
tol = 1e-10;
c = cos(gamma);
s = sin(gamma);

% get axes limits
X = ax.XLim;
Y = ax.YLim;

% create points on line with x limits
if abs(c) > tol
    % we have a non-horizontal line
    P0 = [X; X * s / c];
else
    % we have a horizontal line (assuming through origin)
    P0 = [0, 0; Y];
end

% create points on line with y limits
if abs(s) > tol
    % we have a non-vertical line
    P1 = [Y * c / s; Y];
else
    % we have a vertical line (assuming through origin)
    P1 = [X; 0, 0];
end

% set min limits
if P0(2, 1) >= Y(1)
    X(1) = P0(1, 1);
    Y(1) = P0(2, 1);
else
    X(1) = P1(1, 1);
    Y(1) = P1(2, 1);
end

% set max limits
if P0(2, 2) <= Y(2)
    X(2) = P0(1, 2);
    Y(2) = P0(2, 2);
else
    X(2) = P1(1, 2);
    Y(2) = P1(2, 2);
end
end