function [ h delay ] = extract_h(y,Chans,Fs,Length,Gap,Linear,Over,Fade,Levels,hlen,delay)

    if (nargin<9)  hlen = 0.1;   end;
    if (nargin<10) delay = 0;    end;            % Will estimate it later on 
    
    [Sweep SweepI Offsets] = create_sweep(Fs,Length,Gap,Linear,Over,Fade);
    Gains = 10.^(Levels/20);

    Length = length(Sweep);

    len = Chans * length(Levels) * Length;
    hh = Convolve(SweepI,y);
    hh = hh(min(Offsets(1)-Length/2+delay+(1:len),end),:);          % Put zero time in the middle
    hh = reshape(hh,Length,Chans,length(Levels),[]);                % TIME SPEAKERS GAINS MICS

    if (delay==0)                                           % Calculate earliest signal
        [ peak at ] = max(abs(hilbert(hh(:,:))));
        delay = min(at) + Length/2 - Offsets(1) - Fs/1000;
        hh = hh(delay+1:end,:,:,:);
    end;

    h = zeros(hlen*Fs,Chans,size(y,2),length(Gains),length(Offsets));  % TIME SPEAKERS MICS GAINS HARMONICS
    for (m=1:size(y,2)) 
        for (g=1:length(Gains))
            for (o=1:length(Offsets))
                h(:,:,m,g,o) = Gains(g) * hh(Offsets(o)-Length/2+(1:hlen*Fs),:,g,m);
            end;
        end;
    end; 
    
