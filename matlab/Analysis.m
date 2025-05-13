clear all;

h = extract_file('20250512_SGR_SIGNAL_34CH_SWEEP_48000_5000_1000_100_169_297_LEVELS_20_1.wav');

Speakers = { 'FrontLeft', 'FrontRight', 'FrontCentre', 'Spare1', 'WideLeft', 'WideRight', ...
             'SideHighLeft', 'SideHighRight', 'SideMidLeft','SideMidRight','RearLeft','RearRight',...
             'BackLeft','BackRight','FrontRoofLeft', 'FrontRoofRight', 'RoofLeft', 'RoofRight',...
             'SWL1', 'SWR1', 'SWL2', 'SWR2', 'SWL3', 'SWR3',...
             'SGRLL1', 'SGRLL2', 'SGRLM', 'SGRLH',...
             'SGRRL1', 'SGRRL2', 'SGRRM', 'SGRRH',...
             'KLPL',  'KLPR' };

for (c=1:length(Speakers)) eval([Speakers{c} '=' sprintf('%d',c) ';']); end;


%% Calculate some filters / crossovers etc
addpath('X:\git\matlab');addpath(genpath('X:\git\matlab'))

load('DBX GLENN.cal');
h = h(:,:,1,1);
h = 10*h;

h(:,[SGRLL1 SGRRL1]) = sqrt(10)*h(:,[SGRLL1 SGRRL1]);       % These got turne up 10dB at amp
h(:,[KLPL KLPR]) = 2.5*h(:,[KLPL KLPR]);                    % These were a bit hot (behind counch)


Groups = { { 'ADAM',           [ FrontLeft FrontRight FrontCentre ],   35, 8, 17000, 8 }, ...
           { 'SGR Low',        [ SGRLL1 SGRRL1 ],   35, 8,   800, 8 }, ...
           { 'SGR Mid',        [ SGRLM SGRRM ],  800, 8,  5000, 8 }, ...
           { 'SGR High',       [ SGRLH SGRRH ], 5000, 8, 17000, 8 }, ...
           { 'Sub',            [ SWL1 SWR1 SWL2 SWR2 SWL3 SWR3 ], [],[], 80, 6 },...
           { 'Klipsch',        [ KLPL KLPR ], 35, 8, 17000, 8 }, ...
           { 'JBL',            [ WideLeft : RoofRight ],  80, 6, 17000, 8 } };

LFE   = [ SWL1 SWR1 SWL2 SWR2 SWL3 SWR3 ];
TWEET = [ SGRLH SGRRH ]; 
UNUSED = [ 4 26 30 RoofLeft ];              % NOTE Roof left broken
INV   = [ SGRLM SGRRM ];

%% Create the banded responses and take off the mic response
Fs = 48000;
Fb = [ 2 4 6 8 11:3:47 logspace(log10(50),log10(24000),50)];
H(:,:,:,:,:) = 10*log10(Response(h,48000,0.2,Fb,.002));
H = H - spline(DBX_GLENN(:,1),DBX_GLENN(:,2),Fb)';

G_nom     =  0;         % Nominal level 

%% Calculate trims
H_lev        = mean(H(Fb>500  & Fb<2000,:))      - G_nom;
H_lev(LFE)   = mean(H(Fb>20   & Fb<100,LFE))     - G_nom;
H_lev(TWEET) = mean(H(Fb>2000 & Fb<10000,TWEET)) - G_nom;
h_g        = 10.^(-H_lev/20);
H     = H + repmat(20*log10(h_g),length(Fb),1);

%% Calculate delays
z = abs(hilbert(h));
for (n=1:size(z,2)) delay(n) = sum(cumsum(z(:,n)/max(z(:,n))>0.1)==0); end;
delay = round(max(delay)-delay+1);

for (n=1:size(h,2))
    h_t(:,n) = double((1:max(delay))==delay(n));
end;


%% Create Targets
LPF   = inline('20*log10(1./(1+(F./Fc).^O))','F','Fc','O');
HPF   = inline('20*log10(1./(1+(Fc./(F+.00001)).^O))','F','Fc','O');

L = lines;

for (g=1:length(Groups))
    t = G_nom + LPF(Fb',Groups{g}{5},Groups{g}{6});
    if (~isempty(Groups{g}{3}))
        t = t + HPF(Fb',Groups{g}{3},Groups{g}{4});
    end;
    T(:,Groups{g}{2})=repmat(t,1,length(Groups{g}{2}));
end;

figure(1); clf;
for (g=1:length(Groups))
    semilogx(Fb,T(:,Groups{g}{2}(1)),'color',L(g,:),'LineWidth',2); hold on;
end;
for (g=1:length(Groups))
    semilogx(Fb,mean(H(:,Groups{g}{2}),2),'--','color',L(g,:)); 
end;

axis([10 24000 -40 40]);

%% Create some EQ
Length  = floor((4096-length(h_t)+1)/8)*8/Fs;
TI = T-H;
%TI = TI .* (10.^(T/20) ./ (.1+10.^(T/20)));


%ChannelsToTweak = [FULL MID];
%BandsToTweak    = 1:sum(Fb<100);
%TI(BandsToTweak,ChannelsToTweak) = TI(BandsToTweak,ChannelsToTweak) - (TI(BandsToTweak,ChannelsToTweak)>=3).*(.5*(TI(BandsToTweak,ChannelsToTweak)-3)) - (TI(BandsToTweak,ChannelsToTweak)>=6).*(.25*(TI(BandsToTweak,ChannelsToTweak)-6)) - (TI(BandsToTweak,ChannelsToTweak)>=9).*(.25*(TI(BandsToTweak,ChannelsToTweak)-9));

% Tidy the Sub
ChannelsToTweak = [LFE];
BandsToTweak    = 1:sum(Fb<30);
TI(BandsToTweak,ChannelsToTweak) = TI(BandsToTweak,ChannelsToTweak) - (TI(BandsToTweak,ChannelsToTweak)>=6).*(.5*(TI(BandsToTweak,ChannelsToTweak)-6)) - (TI(BandsToTweak,ChannelsToTweak)>=12).*(.25*(TI(BandsToTweak,ChannelsToTweak)-12)) - (TI(BandsToTweak,ChannelsToTweak)>=18).*(.25*(TI(BandsToTweak,ChannelsToTweak)-18));
TI(:,ChannelsToTweak) = min(6,TI(:,ChannelsToTweak));

% Limit any boost on the SGR woofer
ChannelsToTweak = [ SGRLL1 SGRRL1 ];
TI(:,ChannelsToTweak) = min(TI(:,ChannelsToTweak),0);


% Take off any sharp peaks
t = filter([1 -2 1],1,TI);  
ChannelsToTweak = [1:34];
TI(2:end-1,ChannelsToTweak) = TI(2:end-1,ChannelsToTweak) + (1/2)*(t(3:end,ChannelsToTweak)<3).*t(3:end,ChannelsToTweak);





%%

G = []; 
for (s=1:size(H,2)) 
%    G(:,s) = interp1([0:9 Fb],[TI(1,s)*((0:9)'/10); TI(:,s)],(0:Fs/2),'pchip'); 
    G(:,s) = interp1([0 Fb],[TI(1,s); TI(:,s)],(0:Fs/2),'pchip'); 
end;
G = 10.^([ G; G(end-1:-1:2,:) ]/20);
h_eq = real(ifft(exp(conj(hilbert(log(G))))));
h_eq(:,UNUSED)=0;
h_eq = h_eq(1:Fs*Length,:).*(repmat(  [ones(floor(Fs*3*Length/4),1); cos((0.5:Fs*Length/4)'/Fs/Length*4*pi/2).^2],1,size(H,2)));

h_filt  = Convolve(h_t*diag(h_g), h_eq);
h_filt(:,INV) = - h_filt(:,INV);

% Lets have a look at what we expect


%%
G = interp1(DBX_GLENN(:,1),DBX_GLENN(:,2),0:Fs/2,'pchip')';
G = 10.^([ -G; -G(end-1:-1:2,:) ]/20);
h_miceq = real(ifft(exp(conj(hilbert(log(G))))));
h_miceq = h_miceq(1:Fs*Length,:).*(repmat(  [ones(floor(Fs*7*Length/8),1); cos((0.5:Fs*Length/8)'/Fs/Length*4*pi).^2],1,size(H,2)));

h_out  = Convolve(h(:,:),Convolve(h_filt,h_miceq));
Res = 0.1;
figure(2); clf; 
for (g=1:length(Groups))
    Spectra(h(:,Groups{g}{2}),48000,Res,'color',L(g,:)); hold on; 
end;
axis([10 24000 -40 20]); grid on;

figure(3); clf; 
for (g=1:length(Groups))
    Spectra(h_out(:,Groups{g}{2}),48000,Res,'color',L(g,:)); hold on; 
end;
axis([10 24000 -40 20]); grid on;

figure(4); clf; 
for (g=1:length(Groups))
    Spectra(h_filt(:,Groups{g}{2}),48000,Res,'color',L(g,:)); hold on; 
end;
axis([10 24000 -80 30]); grid on;

%%
figure(5); clf;
Spectra(h_out(:,[FrontLeft FrontRight]),48000,.1,'color',L(1,:)); hold on;
Spectra(h_out(:,[KLPL KLPR]),48000,.1,'color',L(2,:));
Spectra(h_out(:,[SGRLL1 SGRRL1])+h_out(:,[SGRLM SGRRM])+h_out(:,[SGRLH SGRRH]),48000,.1,'color',L(3,:));
axis([10 24000 -40 20]); grid on;

%%

h_deqx = audioread('SGR04112m.wav');
h_deqx = h_deqx(:,[3 5 7]).*[sqrt(10) 1 1]*10.^(-11/20);
Spectra([Convolve(h(:,SGRLL1),h_deqx(:,1)) + Convolve(h(:,SGRLM),h_deqx(:,2)) + Convolve(h(:,SGRLH),h_deqx(:,3)) ],48000,.1,'color',L(4,:));

%%
[b a] = butter(2,800/24000); 
hl = 1.1*impz(conv(b,b),conv(a,a),2048);

[b a] = butter(2,[800 5000]/24000); 
hm = impz(conv(b,b),conv(a,a),2048);

[b a] = butter(2,5000/24000,'high');
hh = .9*impz(conv(b,b),conv(a,a),2048);

Spectra(conv(hl,h(:,SGRLL1))+conv(hm,h(:,SGRLM))+conv(hh,h(:,SGRLH)),48000,.1,'color',L(5,:)); grid on;

%%
Spectra(2*h(:,KLPL),48000,.1,'color',L(6,:)); grid on;



%%
save 20250512_h_filt_050.mat h_filt

