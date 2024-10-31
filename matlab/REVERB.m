clear;
SETUP_20241027

S = [ 1:48 GENELEC ];
M = mod(1:length(S),length(MICS))+1;

rand("seed",0);
M = M(randperm(length(M)));

Rt = [ 0 0.04 * 1.17.^(0:29) ];
Gr = [ 0 10.^((-20:8)/20)];

%%

randn('seed',0);
H0 = randn(200000,length(M));


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
        H = H / sqrt(length(S));
    
        H(:,1:2:48) = Convolve(H(1:end-length(hw)+1,1:2:48),hw);
        H(:,2:2:48) = Convolve(H(1:end-length(ht)+1,1:2:48),ht);
    end;

    
    F={};
    for (s=1:length(S))
        F{end+1} = { M(s), S(s), H(:,s) };
    end;

    file = fopen(sprintf('FILTERS/REVERB/REVERB_%04d.txt',r-1),'wt');
    for (n=1:length(F))
        PrintFilter2(file,n,F{n}{1},F{n}{2},F{n}{3});
        fprintf(file,'\n\n');
    end;
    fclose(file);
end;
