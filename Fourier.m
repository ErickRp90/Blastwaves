function [Datos]=Fourier(Datos)
%Funcion para realizar la transformada de Fourier de las señales de entrada
%y plotear junto con la señal.
%Los valores de salida son los vectores de tiempo y señal cortados. CELDA.
%Las columnas son: 1) ID, 2) sps, 3) tStart, 4) numero de muestras, 
%5) vector de tiempo, 6) vector de señal, 7) vector de frecuencias, 
%8) vector de amplitud, 9) vector de fase.

%Cantidad de señales
Datlen = size(Datos);
Datlen = Datlen(1);

% N = 2^nextpow2(Datlen);

%Se realiza la Tranformada rapida de Fourier (fft).
for i = 1:Datlen
    Datos{i,7} = fft(Datos{i,6})/Datos{i,4};
end

%Espectro de amplitud single-sided.
for i = 1:Datlen
    Datos{i,8} = abs(Datos{i,7});
end

for i = 1:Datlen
    if mod(Datos{i,4}, 2) == 0
        Datos{i,8} = Datos{i,8}(1:(Datos{i,4})/2+1);
    else
        Datos{i,8} = Datos{i,8}(1:(Datos{i,4}-1)/2+1);
    end
end

for i = 1:Datlen
    Datos{i,8}(2:end-1) = 2*Datos{i,8}(2:end-1);
end

% %Filtro Golay para suavizamiento del espectro.
% for i=1:Datlen
%     Datos{i,8}=sgolayfilt(Datos{i,8},3,51);
% end

%Se realiza el eje de frecuencias.
for i=1:Datlen
    Datos{i,7}=Datos{i,2}*(0:(Datos{i,4}/2))/Datos{i,4};
    Datos{i,7}=Datos{i,7}';
end




