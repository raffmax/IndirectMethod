function getAnimationCompassGait(dataIN,gamma,rep,filename)
%getAnimationCompassGait Summary of this function goes here
%   Detailed explanation goes here

figure
grid on; hold on; box on;
axis([-1, 1, -0.5, 1.5])
y_ = [];
foot = [0; 0];
footout = [];
for k = 1:rep
    X  = [dataIN(:,2)+gamma,dataIN(:,2)-dataIN(:,1),...
          dataIN(:,4),dataIN(:,4)-dataIN(:,3)];
    y_ = [y_;X];

    footout = [footout, reshape(repmat(foot,size(X,1),1),2,size(X,1))];
    foot = foot + calc_foot(X(end,:), gamma);

    X  = [dataIN(2:end-1,1)+gamma,dataIN(2:end-1,1)-dataIN(2:end-1,2),...
          dataIN(2:end-1,3),dataIN(2:end-1,3)-dataIN(2:end-1,4)];
    y_ = [y_;X];

    footout = [footout, reshape(repmat(foot,size(X,1),1),2,size(X,1))];
    foot = foot + calc_foot(X(end,:), gamma);
end

% animate
animateCG(y_, footout, gamma, filename);

end

function sw_foot = calc_foot(x, gamma)
L = 1;
theta = x(1);
phi = x(2);
hip = L * [-sin(theta - gamma); cos(theta - gamma)];
sw_foot = hip - L * [sin(phi + gamma - theta); cos(phi + gamma - theta)];
end

