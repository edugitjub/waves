% Carga los sets de eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
num_sujetos=16;
for sujeto=1:num_sujetos
nombre = ['s',int2str(sujeto),'listo.set'];
EEG = pop_loadset('filename',nombre,'filepath','D:\\sanabria\\cleaning\\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
end
eeglab redraw
%% previo
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


%% aves cond

%1:limite_condis(suj,cond)-1
load limite_condis;
c1_ave = zeros(129,N_samples,N_suj);% inicio la matriz
cond=1;
for suj=1:N_suj
    clear c1;
    c1 = ALLEEG(suj,1).data(:,:,1:limite_condis(suj,cond)-1);%entro en sujeto, los tengo todos cargados
    c1mean = mean(c1,3);
    c1_ave(:,:,suj)= c1mean;%
end
%cond2
c2_ave = zeros(129,N_samples,N_suj);% inicio la matriz
cond=2;
for suj=1:N_suj
    clear c2;
    c2 = ALLEEG(suj,1).data(:,:,limite_condis(suj,cond-1):limite_condis(suj,cond)-1);%entro en sujeto, los tengo todos cargados
    c2mean = mean(c2,3);
    c2_ave(:,:,suj)= c2mean;%
end
%cond3
c3_ave = zeros(129,N_samples,N_suj);% inicio la matriz
cond=3;
for suj=1:N_suj
    clear c3;
    c3 = ALLEEG(suj,1).data(:,:,limite_condis(suj,cond-1):limite_condis(suj,cond)-1);%entro en sujeto, los tengo todos cargados
    c3mean = mean(c3,3);
    c3_ave(:,:,suj)= c3mean;%
end
%cond4
c4_ave = zeros(129,N_samples,N_suj);% inicio la matriz
cond=4;
for suj=1:N_suj
    clear c4;
    c4 = ALLEEG(suj,1).data(:,:,limite_condis(suj,cond-1):limite_condis(suj,cond)-1);%entro en sujeto, los tengo todos cargados
    c4mean = mean(c4,3);
    c4_ave(:,:,suj)= c4mean;%
end







%% para comprobar
x=linspace(-1900,896,700);
figure
plot(x, c4_ave(7,:,13))

c1_ave(86,39,16)
b=ALLEEG(6,1).data(7,:,20);%para comprobar
zz = c1_ave(:,:,6);

C1gave = mean(c1_ave,3);
C2gave = mean(c2_ave,3);
C3gave = mean(c3_ave,3);
C4gave = mean(c4_ave,3);

canal = 7;
figure
subplot(2,2,1)
plot(x,C1gave(canal,:),'b')
subplot(2,2,2)
plot(x,C2gave(canal,:),'r')
subplot(2,2,3)
plot(x,C3gave(canal,:),'k')
subplot(224)
plot(x,C4gave(canal,:),'g')


%canal =129;
figure
plot(x,C1gave(canal,:),'b')%400-400
hold on
plot(x,C2gave(canal,:),'r')%400-900
plot(x,C3gave(canal,:),'k')%900-400
plot(x,C4gave(canal,:),'g')%900-900

canal7c214suj = squeeze(c2_ave(7,:,[1 2 4:11 13:16]));%eliminando sujetos 3 y 12
canal7c2mean14suj = mean(canal7c214suj,2);

canal7c414suj = squeeze(c4_ave(7,:,[1 2 4:11 13:16]));
canal7c4mean14suj = mean(canal7c414suj,2);


figure
plot(x,canal7c2mean14suj,'b')
hold on
plot(x,canal7c4mean14suj,'r','LineWidth',.51)

figure

for suj = 1:16
subplot(4,4,suj)
plot(x,canal7c2(:,suj))
hold on
plot(x,canal7c4(:,suj),'r')
plot(x,canal7c1mean,'r','LineWidth',1.5)
set(gca,'ylim',[-10 10])
title(['s' int2str(suj)])
end



%% jhhj

canal = 7;
sujetos_ = 14; %eliminando 3 y 12

c1mean = mean(squeeze(c1_ave(canal,:,[1 2 4:11 13:16])),2);
c3mean = mean(squeeze(c3_ave(canal,:,[1 2 4:11 13:16])),2);
c2mean = mean(squeeze(c2_ave(canal,:,[1 2 4:11 13:16])),2);
c4mean = mean(squeeze(c4_ave(canal,:,[1 2 4:11 13:16])),2);

figure
plot(x,c1mean,'b')%4-4
hold on
plot(x,c2mean,'r','LineWidth',.51)%4-9
plot(x,c3mean,'k')%9-4
plot(x,c4mean,'g')%9-9

%% cambiando linea base
inx = dsearchn(EEG.times',[-1900 -1500]');% Recurso dsearchn , me devuelve la posición de valores de tiempo en en EEG.times, 
%en este caso la linea base

figure
plot (x,c1mean - mean (c1mean(inx(1,1):inx(2,1))))
hold on
plot (x,c2mean - mean (c2mean(inx(1,1):inx(2,1))),'r')
plot (x,c3mean - mean (c2mean(inx(1,1):inx(2,1))),'k')
plot (x,c4mean - mean (c2mean(inx(1,1):inx(2,1))),'g')


%% sujeto 13 que es muy bueno

canal = 129;
c2s13 = c2_ave(canal,:,13);
c4s13 = c4_ave(canal,:,13);

figure
plot(x,c2s13,'b')
hold on
plot(x,c4s13,'r','LineWidth',1)



c1s13 = c1_ave(canal,:,13);
c3s13 = c3_ave(canal,:,13);

figure
plot(x,c1s13,'b')
hold on
plot(x,c3s13,'r','LineWidth',1)


%% figure 13.11
% definitions, selections... LO HAGO POR SUJETOS una figura para cada
% sujeto 

for cond = 1:4
    figure
for sujetoes = 1:16
%cond = 4;
EEG.data = ALLEEG(sujetoes,1).data(:,:,ens_acum(sujetoes,cond):ens_acum(sujetoes,cond+1)-1);
EEG.trials = FIN_ens_ok(sujetoes,cond);

chan2use = 80;%meter un canal de mis datos EDUU
min_freq =  4;
max_freq = 20;
num_frex = 8;

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
eegfft = fft(reshape(EEG.data(chan2use,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);%el reshape: ponerlo todo en una fila

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG. times',[-1800 -1000]');

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

%figure
subplot(4, 4, sujetoes)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-1500 500],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
title('Logarithmic frequency scaling')
end
end

subplot(122)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-1500 500])
title('Linear frequency scaling')
end
%% sf

canal = 129;
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

%% figure 13.11 ADAPTACION EDUARDO MADRID PARA TODOS LOS SUJETOS
% definitions, selections...


%EEG.data = ALLEEG(10,1).data ;
puntos = length (C2_ristra);%
EEG.trials = puntos/700; %ALLEEG(10,1).trials;

chan2use = 80;%meter un canal de mis datos EDUU
min_freq =  4;
max_freq = 40;
num_frex = 30;

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
eegfft = fft(C2_ristra,n_conv_pow2);%el reshape: ponerlo todo en una fila

% initialize
eegpower = zeros(num_frex,EEG.pnts); % frequencies X time X trials
baseidx = dsearchn(EEG.times',[-1800 -1000]');

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
set(gca,'clim',[-3 3],'xlim',[-1900 896],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
title('Logarithmic frequency scaling')

subplot(122)
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-3 3],'xlim',[-1900 896])
title('Linear frequency scaling')
%end



