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

ROT  = [ -180:90:180 ];
ENABLE = logical(zeros(size(MODES,1)*size(ROT,1),max(S)));
for (n=1:size(MODES,1))
    ENABLE(n,ARCH(find(MODES(n,:)==1))) = 1;
    if (sum(ENABLE(n,:))>0)
%        ENABLE(n,SURR) = 1;
        ENABLE(n,[33 37 41 45]) = 1;
        ENABLE(n,ROOF) = 1;    
    end;
end;

if (0)
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


%% CREATING THE SPEAKER REMAP

set = -1;
for (m=1:size(ENABLE,1))
    for (r=ROT)
        set = set+1;
        if (sum(ENABLE(m,:))==0) R=0; 
        else  
            if (r==180 || r==-180) 
                R = Remap([-Xp(1,:); Xp(2,:)]+[r; 0],Yp(:,ENABLE(m,:))); 
            end;
            R = Remap(Xp+[r; 0],Yp(:,ENABLE(m,:)));
        end;
        M = zeros(size(R,1),size(ENABLE,2));
        M(:,ENABLE(m,:)) = R;

        if (~(mean(sum(M'))==1 || mean(sum(M'))==0)) warning ("MISSING SOMETHING FROM INPUT"); keyboard; end;
        
        M(:,SG1_T) = M(:,SG1_W);
        if (sum(sum(M))>0) M(:,SUB) = 1; end;

        file = fopen(sprintf('filters\\ARCH_RESOLUTION\\ARCH_RESOLUTION_%04d.txt',set),'wt');
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
    end;
end;




