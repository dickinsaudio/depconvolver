load 20250512_h_filt_050.mat
h_filt_050 = h_filt;

load h_decorr_trial2.mat
h_d = h_decorr;

h_deqx_sgr_1 = audioread('SGR04112m.wav');
h_deqx = h_deqx_sgr_1(:,[3 5 7]).*[sqrt(10) 1 1]*10.^(-11/20);


%% FOR THE AUTHORING TOOL  4 x 4 = 16 MODES

Speakers = { 'FrontLeft', 'FrontRight', 'FrontCentre', 'Spare1', 'WideLeft', 'WideRight', ...
             'SideHighLeft', 'SideHighRight', 'SideMidLeft','SideMidRight','RearLeft','RearRight',...
             'BackLeft','BackRight','FrontRoofLeft', 'FrontRoofRight', 'RoofLeft', 'RoofRight',...
             'SWL1', 'SWR1', 'SWL2', 'SWR2', 'SWL3', 'SWR3',...
             'SGRLL1', 'SGRLL2', 'SGRLM', 'SGRLH',...
             'SGRRL1', 'SGRRL2', 'SGRRM', 'SGRRH',...
             'KLPL',  'KLPR' };

for (c=1:length(Speakers)) eval([Speakers{c} '=' sprintf('%d',c) ';']); end;

Inputs = {   'FrontLeft', 'FrontRight', 'FrontCentre', 'Spare1', 'WideLeft', 'WideRight', ...
             'SideHighLeft', 'SideHighRight', 'SideMidLeft','SideMidRight','RearLeft','RearRight',...
             'BackLeft','BackRight','FrontRoofLeft', 'FrontRoofRight', 'RoofLeft', 'RoofRight',...
             'Sub1', 'Sub2', 'Sub3', 'Sub4', 'Sub5', 'Sub6',...
             'L1','R1', 'L2', 'R2', 'L3', 'R3', 'L4', 'R4',...
             'FL','FR','C','SUB','LS','RS','LB','RB','LW','RW','LTF','RTF','LTS','RTS','LTB','RTB',...
             'DEQXLS', 'DEQXRS', 'DEQXLL', 'DEQXRL', 'DEQXLM', 'DEQXRM', 'DEQXLH', 'DEQXRH',...
             'I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I7', 'I8', ...
             'MIC01','MIC02','MIC03','MIC04','MIC05','MIC06','MIC07','MIC08',...
             'MIC09','MIC10','MIC11','MIC012','MIC13','MIC14','MIC15','MIC16'...        
         };


for (c=1:length(Inputs)) eval([Inputs{c} '=' sprintf('%d',c) ';']); end;

clear M;
M{length(Speakers),length(Inputs),16}=[];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADAM MODES

M{FrontLeft,L1,1}  = h_filt_050(:,FrontLeft);
M{FrontRight,R1,1} = h_filt_050(:,FrontRight);

M{FrontLeft,L1,4}  = 0.28;
M{FrontRight,R1,4} = 0.28;

M{SWL3,DEQXLS,8} = 1;
M{SWR3,DEQXRS,8} = 1;
M{FrontLeft,DEQXLL,8} = 1;
M{FrontRight,DEQXRL,8} = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SGR MODES

M{SGRLL1,L1, 9} = h_filt_050(:,SGRLL1);
M{SGRRL1,R1, 9} = h_filt_050(:,SGRRL1);
M{SGRLM,L1,  9} = h_filt_050(:,SGRLM);
M{SGRRM,R1,  9} = h_filt_050(:,SGRRM);
M{SGRLH,L1,  9} = h_filt_050(:,SGRLH);
M{SGRRH,R1,  9} = h_filt_050(:,SGRRH);

M{SGRLL1,L1,11} = h_deqx(:,1);
M{SGRRL1,R1,11} = h_deqx(:,1);
M{SGRLM,L1, 11} = h_deqx(:,2);
M{SGRRM,R1, 11} = h_deqx(:,2);
M{SGRLH,L1, 11} = h_deqx(:,3);
M{SGRRH,R1, 11} = h_deqx(:,3);

%%
[b a] = butter(2,800/24000); 
hl = -8*sqrt(.1)*.5*impz(conv(b,b),conv(a,a),2048);

[b a] = butter(2,[800 5000]/24000); 
hm = 1.5*.5*impz(conv(b,b),conv(a,a),2048);

[b a] = butter(2,5000/24000,'high');
hh = .5*impz(conv(b,b),conv(a,a),2048);

%Spectra(conv(hl,h(:,SGRLL))+conv(hm,h(:,SGRLM))+conv(hh,h(:,SGRLH)),48000,.1); grid on;

M{SGRLL1,L1, 12} = hl;
M{SGRRL1,R1, 12} = hl;
M{SGRLM,L1,  12} = hm;
M{SGRRM,R1,  12} = hm;
M{SGRLH,L1,  12} = hh;
M{SGRRH,R1,  12} = hh;

M{SWL3,DEQXLS,16}=1;
M{SWR3,DEQXRS,16}=1;
M{SGRLL1,DEQXLL,16}=1;
M{SGRRL1,DEQXRL,16}=1;
M{SGRLM,DEQXLM,16}=1;
M{SGRRM,DEQXRM,16}=1;
M{SGRLH,DEQXLH,16}=1;
M{SGRRH,DEQXRH,16}=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KLIPSCH MODES

M{KLPL,L1, 17} = h_filt_050(:,KLPL); 
M{KLPR,R1, 17} = h_filt_050(:,KLPR);

M{KLPL,L1, 20} = 2;
M{KLPR,R1, 20} = 2;

M{SWL3, DEQXLS,24}=1;
M{SWR3, DEQXRS,24}=1;
M{KLPL,DEQXLL,24}=1;
M{KLPR,DEQXRL,24}=1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATMOS MIX IN MODES

for (n=1:18) M{n,n,25}=h_filt(:,n); end;      % EQ theatre
for (n=19:24) M{n,n,25}=0.4*h_filt(:,n); end;

M{1,1,26}=h_filt(:,1);                        % ADAM and mix down to Klipch
M{1,3,26}=0.7*h_filt(:,1);
M{2,2,26}=h_filt(:,2);
M{2,3,26}=0.7*h_filt(:,2);
M{1,5,26}=h_filt(:,1);
M{2,6,26}=h_filt(:,2);
M{KLPL,7,26}=h_filt(:,KLPL);
M{KLPL,9,26}=h_filt(:,KLPL);
M{KLPL,11,26}=h_filt(:,KLPL);
M{KLPL,13,26}=h_filt(:,KLPL);
M{KLPL,15,26}=h_filt(:,KLPL);
M{KLPL,17,26}=h_filt(:,KLPL);
M{KLPR,8,26}=h_filt(:,KLPR);
M{KLPR,10,26}=h_filt(:,KLPR);
M{KLPR,12,26}=h_filt(:,KLPR);
M{KLPR,14,26}=h_filt(:,KLPR);
M{KLPR,16,26}=h_filt(:,KLPR);
M{KLPR,18,26}=h_filt(:,KLPR);
M{SWL1,Sub1,26}=.4*h_filt(:,SWL1);
M{SWR1,Sub2,26}=.4*h_filt(:,SWR1);
M{SWL2,Sub3,26}=.4*h_filt(:,SWL2);
M{SWR2,Sub4,26}=.4*h_filt(:,SWR2);
M{SWL3,Sub5,26}=.4*h_filt(:,SWL3);
M{SWR3,Sub6,26}=.4*h_filt(:,SWR3);

M{SGRLL1,1,27}=h_filt(:,SGRLL1);                        % SGR and mix down to Klipch
M{SGRLM, 1,27}=h_filt(:,SGRLM);
M{SGRLH, 1,27}=h_filt(:,SGRLH);
M{SGRRL1,2,27}=h_filt(:,SGRRL1);                        
M{SGRRM, 2,27}=h_filt(:,SGRRM);
M{SGRRH, 2,27}=h_filt(:,SGRRH);
M{SGRLL1,3,27}=0.7*h_filt(:,SGRLL1);                        
M{SGRLM, 3,27}=0.7*h_filt(:,SGRLM);
M{SGRLH, 3,27}=0.7*h_filt(:,SGRLH);
M{SGRRL1,3,27}=0.7*h_filt(:,SGRRL1);                        
M{SGRRM, 3,27}=0.7*h_filt(:,SGRRM);
M{SGRRH, 3,27}=0.7*h_filt(:,SGRRH);
M{SGRLL1,5,27}=h_filt(:,SGRLL1);   
M{SGRLM, 5,27}=h_filt(:,SGRLM);
M{SGRLH, 5,27}=h_filt(:,SGRLH);
M{SGRRL1,6,27}=h_filt(:,SGRRL1);                        
M{SGRRM, 6,27}=h_filt(:,SGRRM);
M{SGRRH, 6,27}=h_filt(:,SGRRH);

M{KLPL,7,27}=h_filt(:,KLPL);
M{KLPL,9,27}=h_filt(:,KLPL);
M{KLPL,11,27}=h_filt(:,KLPL);
M{KLPL,13,27}=h_filt(:,KLPL);
M{KLPL,15,27}=h_filt(:,KLPL);
M{KLPL,17,27}=h_filt(:,KLPL);
M{KLPR,8,27}=h_filt(:,KLPR);
M{KLPR,10,27}=h_filt(:,KLPR);
M{KLPR,12,27}=h_filt(:,KLPR);
M{KLPR,14,27}=h_filt(:,KLPR);
M{KLPR,16,27}=h_filt(:,KLPR);
M{KLPR,18,27}=h_filt(:,KLPR);
M{SWL1,Sub1,27}=.4*h_filt(:,SWL1);
M{SWR1,Sub2,27}=.4*h_filt(:,SWR1);
M{SWL2,Sub3,27}=.4*h_filt(:,SWL2);
M{SWR2,Sub4,27}=.4*h_filt(:,SWR2);
M{SWL3,Sub5,27}=.4*h_filt(:,SWL3);
M{SWR3,Sub6,27}=.4*h_filt(:,SWR3);

for (n=1:18) M{n,n,28}=sqrt(0.1); end;      % Direct theatre
for (n=19:24) M{n,n,28}=0.4*sqrt(0.1); end;     


M{FrontLeft,FL,29}=h_filt(:,FrontLeft);             % Theatre from VID1X1
M{FrontRight,FR,29}=h_filt(:,FrontRight);
M{FrontCentre,C,29}=h_filt(:,FrontCentre);
M{SideMidLeft,LS,29}=h_filt(:,SideMidLeft);
M{SideMidRight,RS,29}=h_filt(:,SideMidRight);
M{BackLeft,LB,29}=h_filt(:,BackLeft);
M{BackRight,RB,29}=h_filt(:,BackRight);
M{WideLeft,LW,29}=h_filt(:,WideLeft);
M{WideRight,RW,29}=h_filt(:,WideRight);
M{FrontRoofLeft,LTF,29}=h_filt(:,FrontRoofLeft);
M{FrontRoofRight,RTF,29}=h_filt(:,FrontRoofRight);
M{SideHighLeft,LTS,29}=h_filt(:,SideHighLeft);
M{SideHighRight,RTS,29}=h_filt(:,SideHighRight);
M{RearLeft,LTB,29}=h_filt(:,RearLeft);
M{RearRight,RTB,29}=h_filt(:,RearRight);
M{SWL1,SUB,29}=.4*h_filt(SWL1);
M{SWL2,SUB,29}=.4*h_filt(SWL2);
M{SWL3,SUB,29}=.4*h_filt(SWL3);
M{SWR1,SUB,29}=.4*h_filt(SWR1);
M{SWR2,SUB,29}=.4*h_filt(SWR2);
M{SWR3,SUB,29}=.4*h_filt(SWR3);


M{FrontLeft,FL,30}=h_filt(:,FrontLeft);             % ADAM + Kippsch from VID
M{FrontRight,FR,30}=h_filt(:,FrontRight);
M{FrontLeft,C,30}=0.7*h_filt(:,FrontLeft);
M{FrontRight,C,30}=0.7*h_filt(:,FrontRight);
M{FrontLeft,LW,30}=h_filt(:,FrontLeft);
M{FrontRight,RW,30}=h_filt(:,FrontRight);
M{KLPL,LS,30}=h_filt(:,KLPL);
M{KLPR,RS,30}=h_filt(:,KLPR);
M{KLPL,LB,30}=h_filt(:,KLPL);
M{KLPR,RB,30}=h_filt(:,KLPR);
M{KLPL,LTF,30}=h_filt(:,KLPL);
M{KLPR,RTF,30}=h_filt(:,KLPR);
M{KLPL,LTS,30}=h_filt(:,KLPL);
M{KLPR,RTS,30}=h_filt(:,KLPR);
M{KLPL,LTB,30}=h_filt(:,KLPL);
M{KLPR,RTB,30}=h_filt(:,KLPR);
M{SWL1,SUB,30}=.4*h_filt(SWL1);
M{SWL2,SUB,30}=.4*h_filt(SWL2);
M{SWL3,SUB,30}=.4*h_filt(SWL3);
M{SWR1,SUB,30}=.4*h_filt(SWR1);
M{SWR2,SUB,30}=.4*h_filt(SWR2);
M{SWR3,SUB,30}=.4*h_filt(SWR3);


M{SGRLL1,FL,31}=h_filt(:,SGRLL1);           % SGR + Kippsch from VID
M{SGRLM, FL,31}=h_filt(:,SGRLM);
M{SGRLH, FL,31}=h_filt(:,SGRLH);
M{SGRRL1,FR,31}=h_filt(:,SGRRL1);                        
M{SGRRM, FR,31}=h_filt(:,SGRRM);
M{SGRRH, FR,31}=h_filt(:,SGRRH);
M{SGRLL1,C ,31}=0.7*h_filt(:,SGRLL1);  
M{SGRLM, C ,31}=0.7*h_filt(:,SGRLM);
M{SGRLH, C ,31}=0.7*h_filt(:,SGRLH);
M{SGRRL1,C ,31}=0.7*h_filt(:,SGRRL1);                        
M{SGRRM, C ,31}=0.7*h_filt(:,SGRRM);
M{SGRRH, C ,31}=0.7*h_filt(:,SGRRH);
M{SGRLL1,LW,31}=h_filt(:,SGRLL1);   
M{SGRLM, LW,31}=h_filt(:,SGRLM);
M{SGRLH, LW,31}=h_filt(:,SGRLH);
M{SGRRL1,RW,31}=h_filt(:,SGRRL1);                        
M{SGRRM, RW,31}=h_filt(:,SGRRM);
M{SGRRH, RW,31}=h_filt(:,SGRRH);

M{KLPL,LS,32}=h_filt(:,KLPL);
M{KLPR,RS,32}=h_filt(:,KLPR);
M{KLPL,LB,31}=h_filt(:,KLPL);
M{KLPR,RB,31}=h_filt(:,KLPR);
M{KLPL,LTF,31}=h_filt(:,KLPL);
M{KLPR,RTF,31}=h_filt(:,KLPR);
M{KLPL,LTS,31}=h_filt(:,KLPL);
M{KLPR,RTS,31}=h_filt(:,KLPR);
M{KLPL,LTB,31}=h_filt(:,KLPL);
M{KLPR,RTB,31}=h_filt(:,KLPR);
M{SWL1,SUB,31}=.4*h_filt(SWL1);
M{SWL2,SUB,31}=.4*h_filt(SWL2);
M{SWL3,SUB,31}=.4*h_filt(SWL3);
M{SWR1,SUB,31}=.4*h_filt(SWR1);
M{SWR2,SUB,31}=.4*h_filt(SWR2);
M{SWR3,SUB,31}=.4*h_filt(SWR3);


M{FrontLeft,FL,32}=sqrt(0.1);                   % Theatre from VID no EQ
M{FrontRight,FR,32}=sqrt(0.1);
M{FrontCentre,C,32}=sqrt(0.1);
M{SideMidLeft,LS,32}=sqrt(0.1);
M{SideMidRight,RS,32}=sqrt(0.1);
M{BackLeft,LB,32}=sqrt(0.1);
M{BackRight,RB,32}=sqrt(0.1);
M{WideLeft,LW,32}=sqrt(0.1);
M{WideRight,RW,32}=sqrt(0.1);
M{FrontRoofLeft,LTF,32}=sqrt(0.1);
M{FrontRoofRight,RTF,32}=sqrt(0.1);
M{SideHighLeft,LTS,32}=sqrt(0.1);
M{SideHighRight,RTS,32}=sqrt(0.1);
M{RearLeft,LTB,32}=sqrt(0.1);
M{RearRight,RTB,32}=sqrt(0.1);
M{SWL1,SUB,32}=0.4*sqrt(0.1);
M{SWL2,SUB,32}=0.4*sqrt(0.1);
M{SWL3,SUB,32}=0.4*sqrt(0.1);
M{SWR1,SUB,32}=0.4*sqrt(0.1);
M{SWR2,SUB,32}=0.4*sqrt(0.1);
M{SWR3,SUB,32}=0.4*sqrt(0.1);











%%

%% CREATE SOME REVERB

load 20250512_Reverb.mat


%% CREATING THE SPEAKER REMAP

ROT=1:(size(h_d,3) + length(Reverb));

set = -1;
for (m=1:size(M,3))
    for (r=ROT)
        set = set+1;
%        file = fopen(sprintf('\\\\10.0.0.192/home/git/depconvolver/matlab/Filters/DEQX_DEMO/DEQX_DEMO_%04d.bin',set),'w');
        file = fopen(sprintf('C:/tmp/DEQX_DEMO_%04d.bin',set),'w');
        fprintf("SET %3d  MODE %2d  ROTATION %3d\n",set,m,r);

        Rev = [];
        if (r>size(h_d,3))
             Rev = Reverb{r-size(h_d,3)};
        end;

        filter=1;
        for (i=1:size(M,2))
            for (o=1:size(M,1))
                if (~isempty(M{o,i,m})) 
                    h = M{o,i,m};
                    if (r<=size(h_d,3))
                        if (i==L1) h = conv(h,h_d(:,1,r)); end;
                        if (i==R1) h = conv(h,h_d(:,2,r)); end;
                    end;

                    fprintf("  INPUT %2d OUTPUT %2d  GAIN %5.3f\n",i,o,20*log10(sum(h.^2)));
                    PrintFilter2(file,filter,i,o,h,true);
                    filter = filter+1;
                end;
            end;
        end;
        if (~isempty(Rev))
            for(f=1:length(Rev))
                fprintf("  INPUT %2d OUTPUT %2d  GAIN %5.3f\n",Rev{f}{1},Rev{f}{2},20*log10(sum(Rev{f}{3}.^2)));
                PrintFilter2(file,filter,Rev{f}{1},Rev{f}{2},Rev{f}{3},true);
                filter = filter+1;
            end;
        end;
                
        fclose(file);
    end;
end;


