%% EQ
clear;
load 20241105_HERVICLE_XOVER.mat


Fs = 48000;


%% INPUTS FROM LISA

ELEVATION = [ zeros(1,56) 30 * ones(1,8)];
AZIMUTH   = [ 0:2.5:25 30:5:50 60:10:300 310:5:330 335:2.5:357.6 0:45:315 ];
LISA      = [ 1:length(ELEVATION) ];
LOW       = find(ELEVATION==0);
HIGH      = find(ELEVATION>0);
MICS      = [ length(LISA)+(1:8) ];

INPUTS = {};
for (n=1:length(LISA)) INPUTS{LISA(n)} = sprintf("LISA E%03d A%03d",ELEVATION(n),floor(AZIMUTH(n))); end;
for (n=1:length(MICS)) INPUTS{MICS(n)} = sprintf("MIC %02d",length(AZIMUTH)+n); end;

%% OUTPUTS

SG1     = [ 1:48 ];                     % Doubled for SG1 - we will duplicate the tweeters later
GENELEC = [ 49 50 51 52 53 54 55  ];
SUB     = [ 56 ];
SG1_W   = SG1(1:2:end);
SG1_T   = SG1(2:2:end);


GENELEV = [ zeros(1,7) ];
GENAZIM = [ -30 0 30 90 135 225 270 ];

SG1ELEV = kron([ zeros(1,16) 0 30 0 30 0 30 0 30 ],[1 1]);
SG1AZIM = kron([ -45 -40 -35 -25 -20 -15 -10 -5 5 10 15 20 25 35 40 45 90 45 135 135 225 225 270 315  ],[1 1]);


ARCH    = [ 1 3 5 49 7 9 11 13 15 50 17 19 21 23 25 51 27 29 31 ];  % Only include every woofer
SURR    = [ 33 37 41 45 52 53 54 55 ];                              % for SG1
ROOF    = [ 35 39 43 47 ];                                          % Duplicate tweeter later

EQ = {};
for (n=1:size(h_filt,2)) EQ(n)={h_filt(:,n)}; end;

OUTPUTS = {};
for (n=1:length(SG1)) 
    OUTPUTS{end+1} = sprintf("SG1 E%03d A%03d LF",SG1ELEV(n),floor(SG1AZIM(n))); 
    OUTPUTS{end+1} = sprintf("SG1 E%03d A%03d HF",SG1ELEV(n),floor(SG1AZIM(n))); 
end;
OUTPUTS(GENELEC) = { "GENELEC LEFT", "GENELEC CENTRE", "GENELEC RIGHT", "GENELEC SURR RIGHT", "GENELEC BACK RIGHT", "GENELEC BACK LEFT", "GENELEC SURR LEFT" };
OUTPUTS(SUB)     = { "SUB 01" };

if (0)
    for (s=1:length(OUTPUTS))
        fprintf("        <txchannel danteId=""%d"" mediaType=""audio""><label>%s</label></txchannel>\n",s,OUTPUTS{s});
    end;
    for (s=1:length(INPUTS))
        fprintf("        <rxchannel danteId=""%d"" mediaType=""audio""><name>%s</name></rxchannel>\n",s,INPUTS{s});
    end;
end; 

%% THE SPEAKERS GEOMETRY FOR CALCULATING PANNINGS
R = 1.0;

Xp  = [ AZIMUTH; ELEVATION ];
X   = R * [  cos(Xp(2,:)/180*pi).*cos(Xp(1,:)/180*pi);
            -cos(Xp(2,:)/180*pi).*sin(Xp(1,:)/180*pi);
             2*sin(Xp(2,:)/180*pi) ];

Xsi = R * [  cos(Xp(2,:)/180*pi).*sin(Xp(1,:)/180*pi);
             cos(Xp(2,:)/180*pi).*cos(Xp(1,:)/180*pi);
             2*sin(Xp(2,:)/180*pi) ];

S   = [ SG1 GENELEC ];
Yp  = [ SG1AZIM GENAZIM; SG1ELEV GENELEV ];
Y   = R * [  cos(Yp(1,:)/180*pi);               % Using Cylinder math for now
             sin(Yp(1,:)/180*pi);
             2*sin(Yp(2,:)/180*pi) ];

figure(1); clf;
plot3(Y(1,SG1_W),Y(2,SG1_W), Y(3,SG1_W),'bo','MarkerSize',15,'MarkerFaceColor','b'); hold on;
plot3(Y(1,SG1_T),Y(2,SG1_T), Y(3,SG1_T),'ro','MarkerSize',5,'MarkerFaceColor','r');
plot3(Y(1,GENELEC),Y(2,GENELEC), Y(3,GENELEC),'gs','MarkerSize',30,'MarkerFaceColor','g');
plot3(X(1,:),X(2,:), X(3,:),'kx','MarkerSize',5,'MarkerFaceColor','k');
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















