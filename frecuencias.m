% Carga los sets de eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
num_sujetos=16;
for sujeto=1:num_sujetos
nombre = ['s',int2str(sujeto),'listo.set'];
EEG = pop_loadset('filename',nombre,'filepath','D:\\sanabria\\cleaning\\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end
eeglab redraw
load('ensayosXcond.mat')
load('FIN_ens_ok.mat')
load('limite_condis.mat')
canal = 7;
ensayo = 2;
sujeto = 1;
muestras = 1:700;
C1 = cell(1,16);
C2 = cell(1,16);
C3 = cell(1,16);
C4 = cell(1,16);
for sujetoes = 1:16 % necesito cells porque hay distinto numero de ensayos por sujeto
C1{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,1:ens_acum(sujetoes,2)-1));
C2{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,ens_acum(sujetoes,2):ens_acum(sujetoes,3)-1));
C3{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,ens_acum(sujetoes,3):ens_acum(sujetoes,4)-1));
C4{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,ens_acum(sujetoes,4):ens_acum(sujetoes,5)-1));
end
N_suj = 16;
N_cond = 4;
N_chan = 129;
N_samples = 700;
%% lk
z=ones(N_suj,1);
ens_onepadded=[z,FIN_ens_ok];
for suj = 1:16
clear a;
a = ens_onepadded(suj,:);
id = cumsum(a);%Indice que me da los valores acumulados de n?mero ensayos para cada suj.
ens_acum(suj,:) =  id;% Guardo valores aqu? de ese indice.
end
limite_condis = ens_acum(:,2:5);%los valores marcan donde empieza la siguiente condición.
%% Para cada condición pongo sujetos muestras ensayos en fila para obtención de potencia x frecuencias
canal = 7;
ensayo = 2;
sujeto = 1;
muestras = 1:700;
C1 = cell(1,16);
C2 = cell(1,16);
C3 = cell(1,16);
C4 = cell(1,16);
for sujetoes = 1:16 % necesito cells porque hay distinto numero de ensayos por sujeto
C1{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,1:ens_acum(sujetoes,2)-1));
C2{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,ens_acum(sujetoes,2):ens_acum(sujetoes,3)-1));
C3{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,ens_acum(sujetoes,3):ens_acum(sujetoes,4)-1));
C4{1,sujetoes} = squeeze(ALLEEG(sujetoes,1).data(canal,muestras,ens_acum(sujetoes,4):ens_acum(sujetoes,5)-1));
end
%lo pongo todo en una fila para correr el fft
for sind = 1:16 % cada una de las 16 matrices en fila
av1 = size(C1{1,sind});
C1{2,sind} = reshape(C1{1,sind},1,av1(1,1) * av1(1,2));
av2 = size(C2{1,sind});
C2{2,sind} = reshape(C2{1,sind},1,av2(1,1) * av2(1,2));
av3 = size(C3{1,sind});
C3{2,sind} = reshape(C3{1,sind},1,av3(1,1) * av3(1,2));
av4 = size(C4{1,sind});
C4{2,sind} = reshape(C4{1,sind},1,av4(1,1) * av4(1,2));
end
C1_ristra = C1{2,1};
for ind = 2:16
C1_ristra = [C1_ristra C1{2,ind}];
end
C2_ristra = C2{2,1};
for ind = 2:16
C2_ristra = [C2_ristra C2{2,ind}];
end
C3_ristra = C3{2,1};
for ind = 2:16
C3_ristra = [C3_ristra C3{2,ind}];
end
C4_ristra = C4{2,1};
for ind = 2:16
C4_ristra = [C4_ristra C4{2,ind}];
end
puntos = length (C4_ristra);%
EEG.trials = puntos/700; %ALLEEG(10,1).trials;
chan2use = 7;%meter un canal de mis datos EDUU
min_freq =  3;
max_freq = 40;
num_frex = 20;
% define wavelet parameters
time = -1.9:1/EEG.srate:.896;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s    =  3./(2*pi*frex); % this line is for figure 13.14
% s    = 10./(2*pi*frex); % this line is for figure 13.14
% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;
% get FFT of data
eegfft = fft(C4_ristra,n_conv_pow2);%el reshape: ponerlo todo en una fila
% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG.times',[-1900 -1600]');
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
figure
subplot(121)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-1200 896],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
title('Logarithmic frequency scaling')
subplot(122)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-1200 896])
title('Linear frequency scaling')
%end
