function drawRABBIT(x, stride, gamma)
    % stride = {1,2}

    trafo = [cos(gamma) -sin(gamma); sin(gamma) cos(gamma)]';

    posShift = -[1,0;0,0]*posHipAUTO(x);

    % colors
    grey = [0.2431 0.2667 0.2980];
    brightblue = [0 0.7451 1];
    darkblue   = [0 0.3176 0.6196];

    % get coordinates of rigid bodies
    posSwingFoot = trafo*(posShift+posSwingFootAUTO(x));
    posSwingKnee = trafo*(posShift+posSwingKneeAUTO(x));
    posStanceFoot = trafo*posShift;
    posStanceKnee = trafo*(posShift+posStanceKneeAUTO(x));
    posHip = trafo*(posShift+posHipAUTO(x));
    posHead = trafo*(posShift+posHeadAUTO(x));

    if stride==1
        % plot stance leg / tibia
        plot([posStanceFoot(1) posStanceKnee(1)], [posStanceFoot(2) posStanceKnee(2)], '-', 'Color', darkblue, 'LineWidth', 6);
        % plot stance leg / femur
        plot([posStanceKnee(1) posHip(1)], [posStanceKnee(2) posHip(2)], '-', 'Color', darkblue, 'LineWidth', 6);
    else
        % plot swing leg / femur
        plot([posHip(1) posSwingKnee(1)], [posHip(2) posSwingKnee(2)], '-', 'Color', darkblue, 'LineWidth', 6);
        % plot swing leg / tibia
        plot([posSwingKnee(1) posSwingFoot(1)], [posSwingKnee(2) posSwingFoot(2)], '-', 'Color', darkblue, 'LineWidth', 6);
    end

    % plot torso
    plot([posHip(1) posHead(1)], [posHip(2) posHead(2)], '-', 'Color', grey, 'LineWidth', 6);

    if stride==1
        % plot swing leg / femur
        plot([posHip(1) posSwingKnee(1)], [posHip(2) posSwingKnee(2)], '-', 'Color', brightblue, 'LineWidth', 6);
        % plot swing leg / tibia
        plot([posSwingKnee(1) posSwingFoot(1)], [posSwingKnee(2) posSwingFoot(2)], '-', 'Color', brightblue, 'LineWidth', 6);
    else
        % plot stance leg / tibia
        plot([posStanceFoot(1) posStanceKnee(1)], [posStanceFoot(2) posStanceKnee(2)], '-', 'Color', brightblue, 'LineWidth', 6);
        % plot stance leg / femur
        plot([posStanceKnee(1) posHip(1)], [posStanceKnee(2) posHip(2)], '-', 'Color', brightblue, 'LineWidth', 6);
    end
end