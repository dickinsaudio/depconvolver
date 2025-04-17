clear all;

MICS = [ 11:26 ];
S = [ 25:36 ];

R = 36;

rand("seed",0);

M(:,1) = mod(0:R-1,length(MICS))+1;
M(:,1) = MICS(M(randperm(R),1));

M(:,2) = mod(0:R-1,length(S))+1;
M(:,2) = S(M(randperm(R),2));


%Rt = [ 0 0.04 * 1.17.^(0:20) ];
%Gr = [ 0 10.^(linespace((-20:0)/20)];
%Rt = Rt(2:2:end);
%Gr = Gr(2:2:end);

Rt = logspace(log10(0.1),log10(1.5),9);
Gr = 10.^(linspace(-15,-3,9)/20);






%%

randn('seed',0);
H0 = randn(200000,length(M));
Fs = 48000;

for (r=1:length(Rt))
    H = H0;
    R_t = Rt(r);
    if (R_t==0) H = zeros(1,length(M));
    else
        a   = log(1e-3)/(Fs*R_t);
        H = H.*(min(.3,exp((1:length(H))*a)))';
        H = H./sqrt(sum(H.^2,1));
        [B A] = butter(4,100/24000,'high');
        H = filter(B,A,H);

        %
        % Air absorption - at x the -6dB freq is FAirAbsorption ~ 37kHz
        % F = FAbs/sqrt(x)
        % This corresponds to a variance of Gaussian of x*log(2)/(2/pi/FA^2)
        % Break that up into four exponential filters
        % clear H; H(142,1)=1; H(round((37/10)^2*142),2)=1; H(round((37/5)^2*142),3)=1; H(end+20,1)=0;

        FAirAbs = 25000;
        SpeedOfSound = 340;
        AirVar = (1:length(H))/Fs*SpeedOfSound*log(2)/(2*pi*pi*(FAirAbs/Fs)^2);
        AirA   = (sqrt(4*(AirVar/4)+1)-1)/2./(AirVar/4);
        s = H(1,:);
        for p=1:2
            for (k=1:length(H))           s=(1-AirA(k))*s + AirA(k)*H(k,:); H(k,:)=s; end;
            s=0; for (k=length(H):-1:1)   s=(1-AirA(k))*s + AirA(k)*H(k,:); H(k,:)=s; end;
        end;

        H = H * Gr(r);
        H = H / sqrt(R);
    end;

    
    F={};
    for (s=1:R)
        F{end+1} = { M(s,1), M(s,2), H(:,s) };
    end;

    Reverb{r}=F;
end;

save 20250415_Reverb Reverb
