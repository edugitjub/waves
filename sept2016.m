
time=(-1900:4:896);
sujetos = 16;
canales = 129;
muestras = 700;

for subject =1: sujetos
CC1(subject,:,:) = mean (C{1,subject},3);
end
 


CC1 = ones(sujetos,canales, muestras);
for suj = 1:16
CC1 = 
end

figure
plot(time,C1m(86,:)
hold on
plot(time,C{1,2}(86,:,1),'r')




