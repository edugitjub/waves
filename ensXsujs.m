% Carga los sets de eeglab
clear all

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

load FIN_ens_ok; load ensayosXcond; 
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

%% sd
C = cell(4,16);%meto en cada fila una de las 4 conds. para cada sujeto canales , muestras (700) y ensayos todos
for cond = 1:4 
    if cond == 1
        for suj=1:N_suj
        C{1,suj}= ALLEEG(suj,1).data(:,:,1:limite_condis(suj,cond)-1);%entro en sujeto, los tengo todos cargados
        end
    else
        for suj=1:N_suj
            C{cond,suj}= ALLEEG(suj,1).data(:,:,limite_condis(suj,cond-1):limite_condis(suj,cond)-1);%entro en sujeto, los tengo todos cargados
        end
    end
end



