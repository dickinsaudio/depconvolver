clear;
Fs     = 48000;
len    = 2048;
N      = 8192;
detail = .05;
sev    = 1;
skew   = .2;
start  = 5;
chans  = 2;

F      = (0:N)/(2*N)*Fs;
bands  = ceil(2 * log2(Fs/2) / detail);
bands  = ceil(bands)

Sev    = linspace(0,1,5);
Sev = 0;

for (s=1:length(Sev))

sev = Sev(s);    
n=47;

rand('seed',n)


A = -bands / lambertw(-1,-bands/exp(1)/N);
b = exp(1)/A; 
x = (0:N);
scale(x<A)  = x(x<A);
scale(x>=A) = A*log(b*x(x>=A));

Fb = exp((0:bands)'*log(Fs/2)/bands);
W  = 1 ./ ( 1 + (start./Fb).^2 );

R = rand(bands,chans) - 0.5;
phase  = [ zeros(1,chans); skew + sev*R];
phase  = W .* phase;
phase  = cumsum(phase);

%phase = phase./std(phase);

clear ph;
for (n=1:size(phase,2))
    ph(:,n) = spline(0:bands,phase(:,n),scale)';                              
end

w = interp1(0:bands,[1 ones(1,bands-1) 0 ],scale,'pchip')';
X = exp(-i*pi/2*ph);
X = fade(X,[0,.01],'hann');
h = ifft([X; conj(X(end-1:-1:2,:))]);
h = h(mod(-len/4+(0:end-1),end)+1,:);
h = h(1:len,:).*tukeywin(len,.2);
[b a] = butter(2,20/24);
h=filter(b,a,h);

if (sev==0) h = 0.*h; h(128,:)=1; end;
h_decorr(:,:,s) = h;


%H = 0*H; H(64,:)=1;

figure(1); clf;
plot(h);
ylim([-.4 .6]);

H = h;
H(10000,1)=0;

figure(2); clf;
subplot(211);

res = .03;
figure(3); clf;
subplot(211);
Spectra([H; zeros(10000,size(H,2))],48000,res,'k'); ylim([-10 10]); hold on;
Spectra((H(:,1)+H(:,2))/2,48000,res,'m');  
Spectra((H(10:end,1)+H(1:end-9,2))/2,48000,res,'b');
Spectra((H(1:end-9,1)+H(10:end,2))/2,48000,res,'r');
ylim([-15,3]);
xlim([100 20000]);
grid on;

subplot(212);
Spectra([H; zeros(10000,size(H,2))],48000,res,'k'); ylim([-10 10]); hold on;
Spectra((H(:,1)+H(:,2))/2,48000,res,'m');  
Spectra((H(12:end,1)+H(1:end-11,2))/2,48000,res,'b');
Spectra((H(1:end-7,1)+H(8:end,2))/2,48000,res,'r');
ylim([-15,3]);
xlim([100 20000]);
grid on;


drawnow;
pause;

end;

h_decorr;
save h_decorr_trial2 h_decorr;

