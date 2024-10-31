
%%  Create a LPF for bass
load 20241020_HERVICLE_XOVER.mat


Fs = 48000;

Fcut = 60;
[b, a] = butter(2,Fcut/(Fs/2),'high');

Fc = 100;
[bl, al] = butter(2, Fc/(Fs/2), 'low');   % Low-pass filter coefficients
[bh, ah] = butter(2, Fc/(Fs/2), 'high');  % High-pass filter coefficients

hl = 10*impz(conv(bl,bl),conv(al,al));
hh = -impz(conv(bh,bh),conv(ah,ah));
hl = conv(hl,impz(b,a));   hl = hl(1:2048);  hl= fade(hl,[0 .1]);

h_8020 = 0.7*h_8020;
hh = conv(hh,impz(b,a));
hw = conv(hh,h_8020(:,1)); hw = hw(1:2048);  hw = fade(hw,[0 .1]);
ht = conv(hh,h_8020(:,2)); ht = ht(1:2048);  ht = fade(ht,[0 .1]); 
hl = fade(hl,[0 .1]);

Spectra([[hl hw ht hl+hw+ht]; zeros(10000,4)]); axis([10 24000 -40 20]); grid on;


%% Create a panning

A = [ -45 -40 -35 -25 -20 -15 -10 -5 5 10 15 20 25 35 40 45 60 90 120 150 210 240 270 300];
a = linspace(-45,300,145);
Mw = zeros(length(A),length(a));
Mt = Mw;
for (n=1:length(a))
    d = A-a(n);
    s = find(d==0);
    if (~isempty(s))
        Mw(s,n)=1;
        Mt(s,n)=1;
    else
        s1 = sum(d<0);
        s2 = s1+1;
        g1 = abs(d(s2));
        g2 = abs(d(s1));
        Mw(s1,n)= g1 / (g2 + g1);
        Mt(s1,n)= Mw(s1,n);
        Mw(s2,n)= g2 / (g2 + g1);
        Mt(s2,n)= Mw(s2,n);
    end;
end;

%%
for (n=1:length(a))
    file = fopen(sprintf('Filters\\PAN_AMP\\PAN_AMP2_%04d.txt',n-1),'wt');
    H = {};
    for (s=1:length(A))
        H{end+1} = { 1, 2*(s-1)+1, 1/length(A)*hl + Mw(s,n)*hw };
        H{end+1} = { 1, 2*(s-1)+2, Mt(s,n)*ht};
    end;
    for (n=1:length(H))
        PrintFilter2(file,n,H{n}{1},H{n}{2},H{n}{3});
        fprintf(file,'\n\n');
    end;
    fclose(file);
end;
