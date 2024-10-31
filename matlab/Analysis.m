clear;

Fs = 48000; 
SweepLength  = 5*Fs;
SweepGap = 1*Fs;
SweepSamples = SweepLength - SweepGap;  
LinearTo = 100;  
OverRun = 169; 
FadeOut = 297;

LinearTime  = floor(((SweepSamples)/(1+log((Fs/2+OverRun)/LinearTo))));               % SweepLength of linear sweep
PhaseDiff   = [linspace(0,LinearTo,LinearTime) logspace(log10(LinearTo),log10(Fs/2+OverRun),SweepSamples-LinearTime)]'/Fs;
Sweep       = sin(2*pi*cumsum(PhaseDiff));
Sweep       = Sweep.*[ones(1,SweepSamples-FadeOut) cos((0.5:FadeOut)/FadeOut/2*pi).^2]';
SweepI      = ifft(1./fft([zeros(50000,1); Sweep; zeros(50000,1)])); 
SweepI      = SweepI(50000+(-1999:SweepSamples+8000));
Offsets     = SweepGap/2 + SweepSamples+2000-round(log([ 1 2:5])*(SweepSamples-LinearTime)/(log(Fs/2+OverRun)-log(LinearTo)));  % Offsets for the non linearities
Sweep       = [ zeros(SweepGap/2,1); Sweep; zeros(SweepGap/2,1) ];


X_direct = audioread('miniDSP_no filter [woofer_tweeter] scaled 1.81E+0  V FSD.wav');
X_harman = audioread('miniDSP_EQ2_HARMAN [woofer_tweeter] scaled 1.81E+0  V FSD.wav');
X_8020   = audioread('miniDSP_EQ2_8020 [woofer_tweeter] scaled 1.81E+0  V FSD.wav');
X_8040   = audioread('miniDSP_EQ2_8040 [woofer_tweeter] scaled 1.81E+0  V FSD.wav');

h_direct = Convolve(X_direct,SweepI);
h_harman = Convolve(X_harman,SweepI);
h_8020   = Convolve(X_8020,  SweepI);
h_8040   = Convolve(X_8040,  SweepI);

T = 218145 + (1:2048);
h_direct = h_direct(T,:);
h_8020   = h_8020(T,:);
h_8020   = fade(h_8020,[0 .25]);

h_8020   = h_8020 / mean(sqrt(sum(h_direct.^2)));

H = {};
for (s=1:24)
    H{end+1} = { s, 2*(s-1)+1, h_8020(:,1)};
    H{end+1} = { s, 2*(s-1)+2, h_8020(:,2)};
end;

file = fopen(sprintf('20241020_HERVICLE_XOVER.txt'),'wt');
for (n=1:length(H))
    PrintFilter2(file,n,H{n}{1},H{n}{2},H{n}{3});
    fprintf(file,'\n\n');
end;
fclose(file);

save 20241020_HERVICLE_XOVER h_8020;



