load 20250415_SGR_h_filt_050.mat
h_filt_050 = h_filt;

load 20250415_SGR_h_filt_150.mat
h_filt_150 = h_filt;

load h_decorr_trial2.mat
h_d = h_decorr;

h_deqx_sgr_1 = audioread('SGR04112m.wav');
h_deqx = h_deqx_sgr_1(:,[3 5 7]).*[10 1 1]*10.^(-13/20);


%% FOR THE AUTHORING TOOL  4 x 4 = 16 MODES

Speakers = { 'ADML','ADMR', ...
             'KLPLL', 'KLPLH', 'KLPRL', 'KLPRH', ...
             'DASL', 'DASR', ...
             'SGRLL', 'SGRLM', 'SGRLH', 'SGRRL', 'SGRRM', 'SGRRH', ...
             'SWL1', 'SWR1', 'SWL2', 'SWR2', 'SWL3', 'SWR3' };

for (c=1:length(Speakers)) eval([Speakers{c} '=' sprintf('%d',c) ';']); end;

Inputs = { 'LEFT', 'RIGHT', ...
           'DEQXLS', 'DEQXRS', 'DEQXLL', 'DEQXRL', 'DEQXLM', 'DEQXRM', 'DEQXLH', 'DEQXRH', ...
         };

for (c=1:length(Inputs)) eval([Inputs{c} '=' sprintf('%d',c) ';']); end;

clear M;
M{length(Speakers),length(Inputs),16}=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADAM MODES

M{ADML,LEFT,1}  = h_filt_050(:,ADML);
M{ADMR,RIGHT,1} = h_filt_050(:,ADMR);

M{ADML,LEFT,2}  = h_filt_150(:,ADML);
M{ADMR,RIGHT,2} = h_filt_150(:,ADMR);

M{ADML,LEFT,7}  = 0.4;
M{ADMR,RIGHT,7} = 0.4;


M{SWL3,DEQXLS,8} = 1;
M{SWR3,DEQXRS,8} = 1;
M{ADML,DEQXLL,8} = 1;
M{ADMR,DEQXRL,8} = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KLIPSCH MODES

M{KLPLL,LEFT,  9} = h_filt_050(:,KLPLL);
M{KLPRL,RIGHT, 9} = h_filt_050(:,KLPRL);
M{KLPLH,LEFT,  9} = h_filt_050(:,KLPLH);
M{KLPRH,RIGHT, 9} = h_filt_050(:,KLPRH);

M{KLPLL,LEFT, 10} = h_filt_150(:,KLPLL);
M{KLPRL,RIGHT,10} = h_filt_150(:,KLPRL);
M{KLPLH,LEFT, 10} = h_filt_150(:,KLPLH);
M{KLPRH,RIGHT,10} = h_filt_150(:,KLPRH);

M{KLPLL,LEFT, 15} = 1;
M{KLPRL,RIGHT,15} = 1;
M{KLPLH,LEFT, 15} = 1;
M{KLPRH,RIGHT,15} = 1;

M{SWL3, DEQXLS,16}=1;
M{SWR3, DEQXRS,16}=1;
M{KLPLL,DEQXLL,16}=1;
M{KLPRL,DEQXRL,16}=1;
M{KLPLH,DEQXLM,16}=1;
M{KLPRH,DEQXRM,16}=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DA MODES

M{DASL,LEFT, 17} = h_filt_050(:,DASL);
M{DASR,RIGHT,17} = h_filt_050(:,DASR);
M{SWL3,LEFT, 17} = h_filt_050(:,SWL3);
M{SWR3,RIGHT,17} = h_filt_050(:,SWR3);

M{DASL,LEFT, 18} = h_filt_150(:,DASL);
M{DASR,RIGHT,18} = h_filt_150(:,DASR);
M{SWL3,LEFT, 18} = h_filt_150(:,SWL3);
M{SWR3,RIGHT,18} = h_filt_150(:,SWR3);


M{DASL, LEFT,  22}=h_filt_050(:,DASL);
M{DASR, RIGHT, 22}=h_filt_050(:,DASR);

M{DASL, LEFT,  23}=-1.8;
M{DASR, RIGHT, 23}=1.8;

M{SWL3,DEQXLS,24}=1;
M{SWR3,DEQXRS,24}=1;
M{DASL,DEQXLL,24}=-1;
M{DASR,DEQXRL,24}=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SGR MODES

M{SGRLL,LEFT,  25} = h_filt_050(:,SGRLL);
M{SGRRL,RIGHT, 25} = h_filt_050(:,SGRRL);
M{SGRLM,LEFT,  25} = h_filt_050(:,SGRLM);
M{SGRRM,RIGHT, 25} = h_filt_050(:,SGRRM);
M{SGRLH,LEFT,  25} = h_filt_050(:,SGRLH);
M{SGRRH,RIGHT, 25} = h_filt_050(:,SGRRH);

M{SGRLL,LEFT,  26} = h_filt_050(:,SGRLL);
M{SGRRL,RIGHT, 26} = h_filt_050(:,SGRRL);
M{SGRLM,LEFT,  26} = h_filt_050(:,SGRLM);
M{SGRRM,RIGHT, 26} = h_filt_050(:,SGRRM);
M{SGRLH,LEFT,  26} = h_filt_050(:,SGRLH);
M{SGRRH,RIGHT, 26} = h_filt_050(:,SGRRH);

M{SGRLL,LEFT,  27} = h_deqx(:,1);
M{SGRRL,RIGHT, 27} = h_deqx(:,1);
M{SGRLM,LEFT,  27} = h_deqx(:,2);
M{SGRRM,RIGHT, 27} = h_deqx(:,2);
M{SGRLH,LEFT,  27} = h_deqx(:,3);
M{SGRRH,RIGHT, 27} = h_deqx(:,3);

%%
[b a] = butter(2,800/24000); 
hl = -8*.5*impz(conv(b,b),conv(a,a),2048);

[b a] = butter(2,[800 5000]/24000); 
hm = 1.5*.5*impz(conv(b,b),conv(a,a),2048);

[b a] = butter(2,5000/24000,'high');
hh = .5*impz(conv(b,b),conv(a,a),2048);

%Spectra(conv(hl,h(:,SGRLL))+conv(hm,h(:,SGRLM))+conv(hh,h(:,SGRLH)),48000,.1); grid on;

M{SGRLL,LEFT,  31} = hl;
M{SGRRL,RIGHT, 31} = hl;
M{SGRLM,LEFT,  31} = hm;
M{SGRRM,RIGHT, 31} = hm;
M{SGRLH,LEFT,  31} = hh;
M{SGRRH,RIGHT, 31} = hh;


%%


M{SWL3,DEQXLS,32}=1;
M{SWR3,DEQXRS,32}=1;
M{SGRLL,DEQXLL,32}=1;
M{SGRRL,DEQXRL,32}=1;
M{SGRLM,DEQXLM,32}=1;
M{SGRRM,DEQXRM,32}=1;
M{SGRLH,DEQXLH,32}=1;
M{SGRRH,DEQXRH,32}=1;

%% CREATE SOME REVERB

load 20250415_Reverb.mat


%% CREATING THE SPEAKER REMAP

ROT=1:(size(h_d,3) + length(Reverb));

set = -1;
for (m=1:size(M,3))
    for (r=ROT)
        set = set+1;
        file = fopen(sprintf('C:\\Tmp\\SGR_DEMO_%04d.bin',set),'w');
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
                        if (i==LEFT)  h = conv(h,h_d(:,1,r)); end;
                        if (i==RIGHT) h = conv(h,h_d(:,2,r)); end;
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


