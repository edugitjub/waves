%OBJETIVO: añadir muestras por los extremos para que no haya artefactos de
%borde. Para datos Sanabria y Correa, expectativa con tonos.
%El input es  cell una de sus filas C 4 filas y 16 columnas conds x sujs 
%Output es una matriz igual que C pero con los extremos reflejados
%para que no haya artefactos de borde.
N_chan = 129;
n_samples = 2100; %añado muestras detras y delante 700 700 700

TCs = cell (4,1); %defino cells
Afl = cell(4,16);

for cond = 1:4
    for suj =1 : 16
        A = C{cond,suj}; % C existe, es una cell con 4 filas y 16 columnas 
        di = size (A);
        n_trials = di (1,3); %en el tercer valor esta la tercera dim de A, ensayos.
        dat_suj = zeros (N_chan,n_samples,n_trials);

        for trial_n = 1:n_trials
            t = squeeze(A (:,:,trial_n));% para cada ensayo...
            T = fliplr (t);% lo giro imagen especular
            ttt = [T t T]; % añado por delante y por detras
            dat_suj(:,:,trial_n) = ttt; % lo meto en la matriz canales x muestras x ensayos
        end
        Afl{cond,suj} = dat_suj(:,576:1525,:); % y para todos los sujetos. Me quedo, de las 700 muestras, con 125 por delante y 125 por detras, o 500 ms adicionales por cada sitio
    end
    % ahora segun la condición. Meto todos los trials para todos los sujetos en
    % la misma dimension. Listo para ir el tfviewerx
    TCs{cond,1} = cat(3, Afl{cond,1}, Afl{cond,2} );
    for sujeto = 3:16 %justo en la linea anterior ya meto los sujetos 1 y 2 así que ahora es a partir del 3
        %cff = cat(3, cff, Afl{cond,sujeto} ); % dim canales x muestras x trials totales en condicion. muestras hay ahora 125 + 700 + 125 = 950
        TCs{cond,1} = cat(3, TCs{cond,1}, Afl{cond,sujeto} );% con 950 muestras
    end
end

