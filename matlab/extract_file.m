function [ h delay Chans Fs Length Gap Linear Over Fade Levels ] = extract_file(file,hlen,delay)

    if (nargin<2) hlen = 0.1;    end;
    if (nargin<3) delay = 0;     end;            % Will estimate it later on 

    Cache = [file(1:end-4) '_cache.mat'];
    if (exist(Cache)) load(Cache); return; end;

    Parts = strsplit(file(1:end-4),'_');
    
    SIGNAL = find(cellfun(@(x) strcmp(x,'SIGNAL'),Parts));  if (isempty(SIGNAL))                         error('Cannot parse file name to extract'); end;
    Chans  = sscanf(Parts{SIGNAL+1},"%d");                  if (isempty(Chans) || Chans<0 || Chans>1024) error('Invalid channels'); end; 

    if (~strcmp(Parts{SIGNAL+2},"SWEEP")) error('Only SWEEP supported at the moment'); end;
    Fs     = sscanf(Parts{SIGNAL+3},"%d");
    Length = sscanf(Parts{SIGNAL+4},"%d");
    Gap    = sscanf(Parts{SIGNAL+5},"%d");
    Linear = sscanf(Parts{SIGNAL+6},"%d");
    Over   = sscanf(Parts{SIGNAL+7},"%d");
    Fade   = sscanf(Parts{SIGNAL+8},"%d");

    if (length(Parts)<SIGNAL+9 || strcmp(Parts{SIGNAL+9},"LEVELS"))
        Levels = 20;
    else
        level  = SIGNAL+10;
        Levels = [];
        while (level<=length(Parts) && ~isempty(sscanf(Parts{level},"%d"))) Levels = [ Levels sscanf(Parts{level},"%d") ]; level=level+1; end;
    end;

    [y fs]      = audioread(file);  if (fs ~= Fs) error('Incorrect sample rate.'); end;
    [ h delay ] = extract_h(y,Chans,Fs,Length,Gap,Linear,Over,Fade,Levels,hlen,delay);

    save(Cache','h','delay','Chans','Fs','Length','Gap','Linear','Over','Fade','Levels','-v6');

