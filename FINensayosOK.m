%INPUTS N_suj N_cond y Ensayos iniciales x condici?n que es matriz N_suj x
%N_cond y ensayos que quito por sujeto que es una cell
% lo que hace este codigo es devolver los ensayos por condici?n una vez
% eliminados los malos FIN_ens_ok. de manera que con eso vuelve uno a los
% datos y hace los aves o lo que haga falta. Enjoy
% Eduardo Madrid, noviembre de 2015
%clear all
n_e = ensayosXcondicion;

%load n_e;%matriz con los ensayos por condici?n, Obtenidos al importar del Net Station, que los voy guardando 
         %puesto que importo una condici?n x file. para cada file muestras
         %totales / muestras x segmento =  ensayos.Hay que a?adir una columna de ceros para que sirve como origen al primer bin 
         %puesto que cuatro intervalos requieren 5 posiciones, entre 5dedos hay 4 espacios  
%load ens_quito_num % deber?a venir de antes, si asi lo hiciera habria que modificar el bucle principal para meter los valores de esta cell 

N_suj = 16;
N_cond = 4;
n_c = N_cond +1;

z=ones(N_suj,1);         
ens_iniXcond_zeropadded = [z,n_e]; % Ensayos iniciales por condici?n.                   
ens_iniXcond_acum = zeros(N_suj,n_c); %son los bins para encajar cuantos de los ensayos eliminados hay en cada bin o condici?n.
num_ens_quitoXcond = zeros (N_suj,N_cond);
%ens_quito_num %aparecen en s1, s2 etc

%% Esto lo puedo hacer al importar de manera que el input a este codigo sea
%% la cell ens_quito_num que tiene tantas filas como sujetos y cada valor
%% es la matriz de los ensayos que quito para ese sujeto.
%ens_quito_num label que quito. Obtenidos al ver eegh justo despues de
%rechazar los ensayos
s1 = [22 33 34 37 48 52 53 54 89 103 108 112 116 117 122 125 148 157 159 176 177 191 192 195 201 202 229 230 243 244 248 260 262 264 276];
s1  = [22 33 34 37 48 52 53 54 89 103 108 112 116 117 122 125 148 157 159 176 177 191 192 195 201 202 229 230 243 244 248 260 262 264 276];
s2  = [44 45 49 54 56 79 106 172 237 244 260];%ensayos rechazados.
s3  = [4 5 34 46 85 102 125 179 189 227 229 240 254 263 264];
s4  = [5 29 66 191 192 235];
s5  = [39 172 178];
s6  = [61 117 136 137 201];
s7  = [5 14 16 36 39 52 57 65 70 91 107 112 133 138 196 203 226 231 233 241 244 250 256 268 275];
s8  = [13 30 38 51 52 69 73 74 115 116 127 131 134 198 205 251 258 277];
s9  = [127 145 148 171 178 186 194 249];
s10  = [34 69 76 99 101 105 106 109 142 171 183 192 195 205 227 235 236];
s11  = [43 48 135 137 138 203 241 266 273 279];
s12  = [6 23 61 66 68 86 111 119 122 136 193 197 202];
s13  = [67 132 155];
s14  = [50 123];
s15  = [10 15 64 96 104 239];
s16  = [0];

ens_quito_num = cell(N_suj,1);% una cell, que ahora no uso, para meter los s1 a s16 pero antes en el proceso.
for suj = 1:N_suj
    matriz=eval(['s',int2str(suj)]);
    ens_quito_num{suj,1}= matriz;
end

%% Bucle principal, el que lo hace.

for suj = 1:16
    clear a;
a = ens_iniXcond_zeropadded(suj,:);
id = cumsum(a);%Indice que me da los valores acumulados de n?mero ensayos para cada suj. 
ens_iniXcond_acum(suj,:) =  id;% Guardo valores aqu? de ese indice. 

matriz=eval(['s',int2str(suj)]);
freq(suj,:) = histc(matriz,id);% saca cuantos ensayos hay para cada bin

end

ens_iniXcond = ens_iniXcond_zeropadded(:,2:5); % le quito la columna de ceros
num_ens_quitoXcond = freq (:,1:4); % idem, no cojo el ultimo bin

FIN_ens_ok = ens_iniXcond - num_ens_quitoXcond;% resta los elementos de dos matrices, GENIAL. 
save FIN_ens_ok FIN_ens_ok
