function [Datos] = ventana2(Datos)
%Funcion para cortar una señal (o vector) para analizarla posteriormente.
%Se necesitan los valores del tiempo inicial y tiempo final.
%Se necesitan, tambien, los valores del vector de tiempo y de la señal para
%cortarlos.CELDA.
%Los valores de salida son los vectores de tiempo y señal cortados. CELDA.
%Las columnas son: 1) ID, 2) sps, 3) tStart, 4) numero de muestras, 
%5) vector de tiempo, 6) vector de señal, 7) vector de frecuencias, 
%8) vector de amplitud, 9) vector de fase.

%Se carga la figura del sismo para seleccionar dos puntos con el mouse para
%realizar el corte de todas las señales.
openfig('Señales.fig');
%ylim([-0.3 0.3])

%% use mouse button to zoom in or out
% Press Enter to get out of the zoom mode
zoom on;

% Wait for the most recent key to become the return/enter key
waitfor(gcf, 'CurrentCharacter', char(13))
zoom reset
zoom off

[x,~] = ginput(1);

%Se busca el valor mas cercano a los datos de entrada del mouse.
Datlen = size(Datos,1);%Numero de filas de la celda. Siempre 3.
for i = 1:Datlen
    [~,idx1] = min(abs(Datos{i,5}-x(1,1)));%Vector con valor del minimo y el indice
    idxin(i,1) = idx1;
end

%Se realiza una nueva celda con los datos del recorte.REAL!!!!!!!!!1!!
for i = 1:Datlen
    Datos{i,4} = 4101;  %Para explosiones 21 s. 1s antes y 20 despues.
    for j = 1:4101
        tven(j,1) = Datos{i,5}(idxin(i,1)-100+j,1);
        sven(j,1) = Datos{i,6}(idxin(i,1)-100+j,1);
    end
    Datos{i,5} = tven;
    Datos{i,6} = sven;
    Datos{i,3} = Datos{i,5}(1,1); %Tiempo de inicio de la señal cortada.
end

