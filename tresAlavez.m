
%% figure 13.11
% definitions, selections...
%chan2use = 7;%meter un canal de mis datos EDUU
clear all
bsl = [-50 0];

load c1; load chanlocs;
EEG.data = c1; EEG.chanlocs = chanlocs;

EEG.srate =250;
EEG.times = -1900:4:896;
EEG.trials = length(EEG.data);

EEG.pnts = 700;

datos =  zeros (129, 15, 700);

for chan2use =1:129
min_freq =  2;
max_freq = 40;
num_frex = 15;

% define wavelet parameters WAVELETT. TENGO DOS ONDAS, LA DE LA SEÑAL Y LA
% WAVELET PARA LA CONVOLUCION
time = -1:1/EEG.srate:1;% En saltos de 1/250 desde el -1 al 1, para la wavelet, centrado en cero.
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;%es la longitud total de los datos mas el vector para la convolución
n_conv_pow2          = pow2(nextpow2(n_convolution));% nextpow2, el siguiente a un núm que es potencia de 2. Ej. para 519 es 10, 2 elevada a 10 es 1024
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(reshape(EEG.data(chan2use,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG.times',bsl');%coloca la linea base

% loop through frequencies and compute synchronization
for fi=1:num_frex
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end

datos (chan2use,:,:) = eegpower;
end
    
tfviewerx(EEG.times, frex, datos, EEG.chanlocs,'c1' );

load gong
soundsc(y)
%% figure 13.11
% definitions, selections...
%chan2use = 7;%meter un canal de mis datos EDUU
clear all
load c2; load chanlocs;
EEG.data = c2; EEG.chanlocs = chanlocs;

EEG.srate =250;
EEG.times = -1900:4:896;
EEG.trials = length(EEG.data);

EEG.pnts = 700;

datos =  zeros (129, 15, 700);

for chan2use =1:129
min_freq =  2;
max_freq = 40;
num_frex = 15;

% define wavelet parameters WAVELETT. TENGO DOS ONDAS, LA DE LA SEÑAL Y LA
% WAVELET PARA LA CONVOLUCION
time = -1:1/EEG.srate:1;% En saltos de 1/250 desde el -1 al 1, para la wavelet, centrado en cero.
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;%es la longitud total de los datos mas el vector para la convolución
n_conv_pow2          = pow2(nextpow2(n_convolution));% nextpow2, el siguiente a un núm que es potencia de 2. Ej. para 519 es 10, 2 elevada a 10 es 1024
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(reshape(EEG.data(chan2use,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG.times',[0 100]');%coloca la linea base

% loop through frequencies and compute synchronization
for fi=1:num_frex
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end

datos (chan2use,:,:) = eegpower;
end

tfviewerx(EEG.times, frex, datos, EEG.chanlocs,'c2' );


%% figure 13.11
% definitions, selections...
%chan2use = 7;%meter un canal de mis datos EDUU
clear all
load c3; load chanlocs;
EEG.data = c3; EEG.chanlocs = chanlocs;

EEG.srate =250;
EEG.times = -1900:4:896;
EEG.trials = length(EEG.data);

EEG.pnts = 700;

datos =  zeros (129, 15, 700);

for chan2use =1:129
min_freq =  2;
max_freq = 40;
num_frex = 15;

% define wavelet parameters WAVELETT. TENGO DOS ONDAS, LA DE LA SEÑAL Y LA
% WAVELET PARA LA CONVOLUCION
time = -1:1/EEG.srate:1;% En saltos de 1/250 desde el -1 al 1, para la wavelet, centrado en cero.
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;%es la longitud total de los datos mas el vector para la convolución
n_conv_pow2          = pow2(nextpow2(n_convolution));% nextpow2, el siguiente a un núm que es potencia de 2. Ej. para 519 es 10, 2 elevada a 10 es 1024
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(reshape(EEG.data(chan2use,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG.times',[0 100]');%coloca la linea base

% loop through frequencies and compute synchronization
for fi=1:num_frex
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end

datos (chan2use,:,:) = eegpower;
end


    
tfviewerx(EEG.times, frex, datos, EEG.chanlocs,'c3' );




%% figure 13.11
% definitions, selections...
%chan2use = 7;%meter un canal de mis datos EDUU
clear all
load c4; load chanlocs;
EEG.data = c4; EEG.chanlocs = chanlocs;

EEG.srate =250;
EEG.times = -1900:4:896;
EEG.trials = length(EEG.data);

EEG.pnts = 700;

datos =  zeros (129, 15, 700);

for chan2use =1:129
min_freq =  2;
max_freq = 40;
num_frex = 15;

% define wavelet parameters WAVELETT. TENGO DOS ONDAS, LA DE LA SEÑAL Y LA
% WAVELET PARA LA CONVOLUCION
time = -1:1/EEG.srate:1;% En saltos de 1/250 desde el -1 al 1, para la wavelet, centrado en cero.
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;%es la longitud total de los datos mas el vector para la convolución
n_conv_pow2          = pow2(nextpow2(n_convolution));% nextpow2, el siguiente a un núm que es potencia de 2. Ej. para 519 es 10, 2 elevada a 10 es 1024
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(reshape(EEG.data(chan2use,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG.times',[0 100]');%coloca la linea base

% loop through frequencies and compute synchronization
for fi=1:num_frex
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end

datos (chan2use,:,:) = eegpower;
end


    
tfviewerx(EEG.times, frex, datos, EEG.chanlocs,'c4' );


load gong
soundsc(y)


