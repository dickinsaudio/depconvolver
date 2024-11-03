SETUP_20241027;

%% FOR THE AUTHORING TOOL  4 x 4 = 16 MODES


ZOOM = [0.5 1.0 2.0 3.0];
%ZOOM = [ 1.0 ];
ROT  = [ -180:5:180 ];

EN_27_4 = [ ARCH(logical([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1])) SURR ROOF ];
EN_23_4 = [ ARCH(logical([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1])) 33 37 41 45 ROOF ];
EN_7_4  = [ GENELEC ROOF ];
EN_7_4  = [ GENELEC(1:3) 33 37 41 45 ROOF ];

ENABLE = logical(zeros(1,max(S)));
ENABLE(1,EN_23_4) = 1;                  % Full set of speakers
ENABLE(2,:)       = ENABLE(1,:);        % This will be mirrored
ENABLE(3,:)       = 0;                  % Nothing
ENABLE(4,EN_7_4)  = 1;                  % A comparative downmix
ENABLE = kron(ENABLE,ones(size(ZOOM,2),1));
ENABLE = logical(ENABLE);

if (0)
    for (n=1:size(ENABLE,1))
        figure(2); clf;
        plot3(Y(1,ENABLE(n,:)),Y(2,ENABLE(n,:)), Y(3,ENABLE(n,:)),'bo','MarkerSize',15,'MarkerFaceColor','b'); hold on;
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

%% CREATING THE SPEAKER REMAP

set = -1;
for (m=1:size(ENABLE,1))
    for (r=ROT)
        set = set+1;
        if (sum(ENABLE(m,:))==0) R=0; 
        else  
            if (m>=5 && m<=8)
                R = Remap([180-Xp(1,:); Xp(2,:)]+[r; 0],Yp(:,ENABLE(m,:))); 
            else
                R = Remap(Xp+[r; 0],Yp(:,ENABLE(m,:)));
            end;
        end;
        M = zeros(size(R,1),size(ENABLE,2));
        M(:,ENABLE(m,:)) = R;

        if (mean(sum(M'))~=1) keyboard; end;
        
        M(:,SG1_T) = sqrt(M(:,SG1_W));
        M(:,SUB) = 1;

        file = fopen(sprintf('Filters\\AUTHORING\\AUTHORING_%04d.txt',set),'wt');
        fprintf("SET %3d  MODE %2d  ROTATION %3d\n",set,m,r);
        filter=1;
        for (i=1:size(M,1))
            for (o=1:size(M,2))
                if (M(i,o)~=0) 
                    % fprintf("  INPUT %2d OUTPUT %2d  GAIN %5.3f\n",i,o,M(i,o));
                    PrintFilter2(file,filter,i,o,M(i,o)*EQ{o});
                    fprintf(file,'\n\n');
                    filter = filter+1;
                end;
            end;
        end;
        fclose(file);
        set = set + 1;
    end;
end;



