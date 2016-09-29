 %% adaptado de figure 13.11. INPUT es TCs, cell con 4 filas, cada una una condicion con 129 x 950 x trials
EEG.srate =250; EEG.times = -2400:4:1396; min_freq =  2; max_freq = 40; num_frex = 15;EEG.pnts = 950;load chanlocs;EEG.chanlocs = chanlocs;
dataEEGpower = cell (4,1);
dataEEGpowerSubject = cell (4,16);

%for subj =1:16 % por sujetos , quitar este bucle si se hace gave
for cond = 1:4
    %EEG.data = Afl{cond,subj}; %uncom para suj    por sujetos, para estadistica
   EEG.data = TCs{cond,1}; %TCs debe existir, es una cell de dim 4 x 1 para total. 
   
tr = size(EEG.data); EEG.trials = tr(1,3);% numero de ensayos para una cond para un sujeto.
datos =  zeros (129, 15, 950);
      for chan2use =1:129
        % define wavelet parameters WAVELETT. TENGO DOS ONDAS: SEÑAL Y WAVELET PARA LA CONVOLUCION
        time = -1:1/EEG.srate:1;% En saltos de 1/250 desde el -1 al 1, para la wavelet, centrado en cero.
        frex = logspace(log10(min_freq),log10(max_freq),num_frex);
        s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
 
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
      
%dataEEGpowerSubject{cond,subj} = datos;% umcomment this line for sujeto a sujeto para estadistica     
dataEEGpower{cond,1} = datos; %umcomment this line for todos

end %fin final de las 4 cond para sujeto subj

%end % uncomment para hacerlo x sujetos

load gong
soundsc(y)


%% df

cond = 2;
for s =1:16
Subject = s;
datos = dataEEGpowerSubject{cond,Subject}; 
tfviewerx(EEG.times, frex, datos, EEG.chanlocs,['c' int2str(cond) 's' int2str(Subject) ] );
end


for cond = 1:4
datos = dataEEGpower{cond,1};
tfviewerx(EEG.times, frex, datos, EEG.chanlocs,['c' int2str(cond)] );
end

figure
plot (squeeze(dataEEGpower{2,1}(11,5,:)))
%% medias por sujeto. Se necesita la cell con las matrices por sujeto.
for si = 1:16
    for ci = 1:4
v1m = dataEEGpowerSubject{ci,si}([106 113 81 129], 4:6, 480:540);
v1(ci,si) = mean(mean(mean (v1m,3),2));
    end
end

