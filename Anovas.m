%% ANOVAS

%% Datos preparando. medias por sujeto. Se necesita la cell con las matrices por sujeto.
canales = [106 113 81 129];
frecuencias = 4:6;
muestras = 480:540;

for si = 1:16
    for ci = 1:4
v1m = dataEEGpowerSubject{ci,si}(canales, frecuencias, muestras);
v1(ci,si) = mean(mean(mean (v1m,3),2));
    end
end


%introducimos la codificacion del dise?o experimental
rr(1:32,1)=1;
rr(33:64,1)=2;

ict(1:16,1)=1;
ict(17:32,1)=2;
ict(33:48,1)=1;
ict(49:64,1)=2;

sujetos(1:16,1)=1:16;
sujetos(17:32,1)=1:16;
sujetos(33:48,1)=1:16;
sujetos(49:64,1)=1:16;

nombres={'rr', 'ict'};% los nombres de las variables

   
   %datos(1:16,1)= v1(1,:);
   %datos(17:32,1)= v1(2,:);
   %datos(33:48,1)= v1(3,:);
   %datos(49:64,1)= v1(4,:); % la siguiente linea hace lo mismo que estas
   %4. Ojo con las dim de las matrices.
   
   datos(:,1) = reshape(v1',1,64);
      
 stats = rm_anova2(datos,sujetos,rr,ict,nombres);%hace el ANOVA
     
    signifs(canal,1) = canal(1,1);%canal     
    signifs(canal,2) = stats{2,6};%tarea
    signifs(canal,3) = stats{3,6};%cambio
    signifs(canal,4) = stats{4,6};%interaccion
      
