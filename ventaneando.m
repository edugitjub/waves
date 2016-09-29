
%% Cohen code para sacar ventanas de datos


tidx = dsearchn(EEG.times',TimeRange'); % indices, busca las posiciones en EEG.times en las que estan los valores TimeRange.
fidx = dsearchn(frequencies',[freqL freqH]'); % frecuencias. Ojo con lo de traspuestas.

TF_Range_Power = mean(mean( tf_data(fidx(1):fidx(2),tidx(1):tidx(2)) ,1),2); %media de sub-matriz o ventana de tiempos-frecuencias



TimeRange = [2000 2500];
EEG.times = 0:4:3000;