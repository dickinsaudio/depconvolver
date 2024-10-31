SETUP_20241027;

%% FOR THE AURICLE INVESTIGATION - 4 x 8 = 32 MODES

MODES   = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;    % 19          % 45 Degree Width
            1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1;    % 10
            1 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 0 0 1;    %  5
            1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1;    %  3
            1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;    %  2
            0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;    %  1

            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;

            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0;    % 13          % 20 Degree Width
            0 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 0 0;    %  7
            0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0;    %  5
            0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0;    %  3;
            0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;    %  2
            0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;    %  1

            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;     
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;

            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;    % 19          % Shrink around 3 columns
            0 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 1 1 0;    % 15
            0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0;    %  9
            0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0;    %  3
            0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;    %  2
            0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;    %  1

            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;


            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;    % 19          % Narrowing in around centre
            0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0;    % 17
            0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0;    % 15
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0;    % 13
            0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0;    % 11
            0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0;    % 09
            0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0;    % 07
            0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0; ]; % 05

ENABLE = logical(zeros(size(MODES,1),max(S)));
for (n=1:size(MODES,1))
    ENABLE(n,ARCH(find(MODES(n,:)==1))) = 1;
    ENABLE(n,SURR) = 1;
    ENABLE(n,ROOF) = 1;    
end;

if (1)
    for (n=1:size(ENABLE,1))
        figure(2); clf;
        plot3(Y(1,find(ENABLE(n,:)==1)),Y(2,find(ENABLE(n,:)==1)), Y(3,find(ENABLE(n,:)==1)),'bo','MarkerSize',15,'MarkerFaceColor','b'); hold on;
        ax = gca;
        ax.PlotBoxAspectRatio = [1 1 1];
        ax.DataAspectRatioMode = 'manual';
        ax.PlotBoxAspectRatioMode = 'manual';
        ax.CameraPositionMode = 'manual';
        ax.CameraTargetMode = 'manual';
        ax.CameraViewAngleMode = 'manual';
        ax = gca;
        ax.CameraPosition = [30 -10 8];
        ax.CameraTarget = [0 0 0];
        grid on;
        pause;
    end;
end;
