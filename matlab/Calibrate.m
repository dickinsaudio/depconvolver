clear all;

Signal;

X = audioread('20241101_COMO_SETUP_SUB.wav');
X = reshape(X,SweepLength,[]);

Delay = 1200;

h = Convolve(X,SweepI);
h = h(Offsets(1)+Delay+(1:4800),1:56);
h(48000,1)=0;

%%
FULL      =  49:51;
LFE       =  56;
MID       =  1:2:48;
HIGH      =  2:2:48;
OFF       =  52:55;

Fs = 48000;
Fb = [ 5 10 15 logspace(log10(20),log10(24000),50)];
H = 10*log10(Response(h,48000,1,Fb,0));

GAIN      = median(median(H(30:40,:)));

G_nom     =  0;         % Nominal level 
G_lfe     =  0;

H_lev       = median(H(Fb>500   & Fb<2000,:  ))  - G_nom;
H_lev(LFE)  = median(H(Fb>20    & Fb<100 ,LFE))  - G_lfe;
H_lev(HIGH) = median(H(Fb>1000  & Fb<4000,HIGH)) - G_nom; 
h_g         = 10.^(-H_lev/20);
h_g(OFF)=0;
H     = H + repmat(20*log10(h_g),length(Fb),1);
h_g   = h_g * 10.^(GAIN/20);

%% Calculate delays
z = abs(hilbert(h));
for (n=1:size(z,2)) delay(n) = sum(cumsum(z(:,n)/max(z(:,n))>0.1)==0); end;
delay = round(max(delay)-delay+1);
delay(1:2:48) = delay(2:2:48);
delay(OFF) = 0;

for (n=1:size(h,2))
    h_t(:,n) = double((1:max(delay([FULL MID LFE])))==delay(n));
end;

%% Create Targets
LPF   = inline('20*log10(1./(1+(F./Fc).^O))','F','Fc','O');
HPF   = inline('20*log10(1./(1+(Fc./(F+.00001)).^O))','F','Fc','O');
F_top  = 18000;     O_top   = 4;
F_lfe  = 110;       O_lfe   = 3;
F_full = 110;       O_full  = 3; 
F_mid  = 110;       O_mid   = 3;
F_high = 900;       O_high  = 5;

T         =   G_nom + repmat(LPF(Fb',F_top, O_top ),1,size(H,2));
T(:,LFE)  =   G_lfe + repmat(LPF(Fb',F_lfe, O_lfe ),1,length(LFE));
T(:,LFE)  =   T(:,LFE) + repmat(HPF(Fb',15,6 ),1,length(LFE));
T(:,FULL) =   T(:,FULL)+ repmat(HPF(Fb',F_full,O_full),1,length(FULL));
T(:,MID ) =   T(:,MID )+ repmat(HPF(Fb',F_mid ,O_mid ),1,length(MID));
T(:,MID ) =   T(:,MID )+ repmat(LPF(Fb',10000,4),1,length(MID));
T(:,HIGH ) =   T(:,HIGH)+ repmat(HPF(Fb',F_high ,O_high ),1,length(HIGH));


figure(1); clf;
semilogx(Fb,mean(H(:,FULL),2),'b'); hold on;
semilogx(Fb,mean(H(:,MID),2),'m');
semilogx(Fb,mean(H(:,HIGH),2),'r');
semilogx(Fb,mean(H(:,LFE),2),'k');
semilogx(Fb,mean(T(:,FULL),2),'b--');
semilogx(Fb,mean(T(:,HIGH),2),'r--');
semilogx(Fb,mean(T(:,MID),2),'m--');
semilogx(Fb,mean(T(:,LFE),2),'k--');
axis([10 24000 -40 40]);

%% Create some EQ
TI = T-H;
TI(:,OFF)=0;
TI = TI .* (10.^(max(T,H)/20) ./ (.1+10.^(max(T,H)/20)));

%%
ChannelsToTweak = [FULL MID];
BandsToTweak    = 1:sum(Fb<100);
%TI(BandsToTweak,ChannelsToTweak) = TI(BandsToTweak,ChannelsToTweak) - (TI(BandsToTweak,ChannelsToTweak)>=3).*(.5*(TI(BandsToTweak,ChannelsToTweak)-3)) - (TI(BandsToTweak,ChannelsToTweak)>=6).*(.25*(TI(BandsToTweak,ChannelsToTweak)-6)) - (TI(BandsToTweak,ChannelsToTweak)>=9).*(.25*(TI(BandsToTweak,ChannelsToTweak)-9));

ChannelsToTweak = [LFE];
BandsToTweak    = 1:sum(Fb<30);
TI(BandsToTweak,ChannelsToTweak) = TI(BandsToTweak,ChannelsToTweak) - (TI(BandsToTweak,ChannelsToTweak)>=6).*(.5*(TI(BandsToTweak,ChannelsToTweak)-6)) - (TI(BandsToTweak,ChannelsToTweak)>=12).*(.25*(TI(BandsToTweak,ChannelsToTweak)-12)) - (TI(BandsToTweak,ChannelsToTweak)>=18).*(.25*(TI(BandsToTweak,ChannelsToTweak)-18));
TI(:,ChannelsToTweak) = min(10,TI(:,ChannelsToTweak));

t = filter([1 -2 1],1,TI);  
ChannelsToTweak = [1:56];
%TI(2:end-1,ChannelsToTweak) = TI(2:end-1,ChannelsToTweak) + (1/2)*(t(3:end,ChannelsToTweak)<3).*t(3:end,ChannelsToTweak);

TI(end,:)=TI(end-1,:);

%%

G = []; 
for (s=1:size(H,2)) 
%    G(:,s) = interp1([0:9 Fb],[TI(1,s)*((0:9)'/10); TI(:,s)],(0:Fs/2),'pchip'); 
    G(:,s) = interp1([0 Fb],[0; TI(:,s)],(0:Fs/2),'pchip'); 
end;
G = 10.^([ G; G(end-1:-1:2,:) ]/20);
h_eq = real(ifft(exp(conj(hilbert(log(G))))));
h_filt = Convolve(h_t*diag(h_g), h_eq);
h_filt = fade(h_filt(1:2048,:),[0 128]);


%%
%%
h_out  = Convolve(h(:,:),h_filt);
Res = 0.1;
figure(1); 
clf; 
Spectra(h(:,FULL),48000,Res,'b'); hold on; 
Spectra(h(:,MID),48000,Res,'g'); 
Spectra(h(:,HIGH),48000,Res,'r'); 
Spectra(h(:,LFE),48000,Res,'m'); 
axis([10 24000 -80 0]); 
grid on;


figure(2); 
clf; 
Spectra(h_out(:,FULL),48000,Res,'b'); hold on; 
Spectra(h_out(:,MID),48000,Res,'g'); 
Spectra(h_out(:,HIGH),48000,Res,'r'); 
Spectra(h_out(:,LFE),48000,Res,'m'); 
plot(Fb, T+GAIN,'k:');
axis([10 24000 -80 0]); grid on;

figure(3); 
clf; 
Spectra(h_filt(:,FULL),48000,'b'); hold on; 
Spectra(h_filt(:,MID),48000,'g'); 
Spectra(h_filt(:,HIGH),48000,'r'); 
Spectra(h_filt(:,LFE),48000,'m'); 
axis([10 24000 -20 20]); grid on;

%% Now to add a crossover
Fc = 3000;
[bl, al] = butter(1, .7*Fc/(Fs/2), 'low');   % Low-pass filter coefficients
[bh, ah] = butter(1, 1.4*Fc/(Fs/2), 'high');  % High-pass filter coefficients

%hl = impz(conv(bl,bl),conv(al,al));
%hh = impz(conv(bh,bh),conv(ah,ah));
hl = impz(bl,al,256);
hh = -impz(bh,ah,256);

h_filt(:,1:2:48) = filter(hl,1,h_filt(:,1:2:48));
h_filt(:,2:2:48) = filter(hh,1,h_filt(:,2:2:48));

%% Now repeat with the 2 way speakers
h_out = Convolve(h(:,:),h_filt);
h2 = h_out(:,1:2:48) + h_out(:,2:2:48);

H = 10*log10(Response(h2,48000,1,Fb,0));
H_lev       = median(H(Fb>500   & Fb<2000,:  ))  - G_nom;
h_g         = 10.^(-H_lev/20);
H     = H + repmat(20*log10(h_g),length(Fb),1);
h_g   = h_g * 10.^(GAIN/20);

TI = T(:,49)-H;
TI(:,OFF)=0;
TI(end,:)=TI(end-1,:);
TI = TI .* (10.^(T(:,49)/20) ./ (.1+10.^(T(:,49)/20)));
TI = filtfilt([.5 .5],1,TI);

G = []; 
for (s=1:size(H,2)) 
%    G(:,s) = interp1([0:9 Fb],[TI(1,s)*((0:9)'/10); TI(:,s)],(0:Fs/2),'pchip'); 
    G(:,s) = interp1([0 Fb],[0; TI(:,s)],(0:Fs/2),'pchip'); 
end;
G = 10.^([ G; G(end-1:-1:2,:) ]/20);
h_eq = real(ifft(exp(conj(hilbert(log(G))))));
h_filt2 = h_g.*h_eq;
h_filt2 = fade(h_filt2(1:2048,:),[0 128]);

for (n=2:2:48)
    h_filt(:,n-1)=filter(h_filt2(:,n/2),1,h_filt(:,n-1));
    h_filt(:,n  )=filter(h_filt2(:,n/2),1,h_filt(:,n  ));
end;

h_filt = fade(h_filt,[0 256]);


%%

h_out  = Convolve(h(:,:),h_filt);

figure(4);
clf;
Spectra(h_out(:,1)+h_out(:,2),48000,Res,'g'); hold on;
Spectra(h_out(:,1),48000,Res,'b'); hold on;
Spectra(h_out(:,2),48000,Res,'r'); hold on;
axis([10 24000 -80 00]); grid on;

%%
figure(5);
clf;
Spectra(h_out(:,50),48000,Res,'k'); hold on;
Spectra(h_out(:,15)+h_out(:,16),48000,Res,'b'); hold on;
Spectra(h_out(:,17)+h_out(:,18),48000,Res,'r'); hold on;
Spectra(h_out(:,15)+h_out(:,16)+h_out(:,17)+h_out(:,18),48000,Res,'m'); hold on;
Spectra(h_out(:,51)+h_out(:,29)+h_out(:,30),48000,Res,'g'); hold on;
Spectra(h_out(:,51)+h_out(:,49),48000,Res,'c'); hold on;
axis([10 24000 -80 20]); grid on;

%%
save 20241105_HERVICLE_XOVER.mat h_filt





