%Programa para procesar las explosiones registradas en el proyecto
%Huelitli-3D.
clear  clc

%% 
%Leer archivos las tres señales en formato gcf. Se cargan con el
%comando uigetfile en una celda. Se seleccionan todos los archivos.
%Archivos de 1 hora para que no se sobreescriban.
%Falla cuando cargo solo un archivo.

S=uigetfile('*.gcf','Selecciona los archivos','Multiselect','on');%uigetfile-Abre un archivo con cuadro de dialogo.
a=length(S); %Siempre igual a 3.

%La información de las señales se guarda por columnas en una celda.
%El numero de filas son las señales para procesar.
%Las columnas son: 1) ID, 2) sps, 3) tStart, 4) numero de muestras, 
%5) vector de tiempo, 6) vector de señal, 7) vector de frecuencias, 
%8) vector de amplitud, 9) vector de fase.
INFO=cell(a,9);

%Se leen los datos de los archivos cargados con uigetfile.
Data=cell(a,1);%Guardar los datos en una celda.
for i=1:a
    Si=cell2mat(S(i));%Convert cell array to ordinary array of the underlying data type.
    [samples,streamID,sps,tStart]=readgcffile(num2str(Si));
    INFO{i,1}=streamID;
    INFO{i,2}=sps;
    INFO{i,3}=tStart;
    INFO{i,4}=length(samples);
    INFO{i,6}=samples;
end

%%
% Pre-pocesamiento de la informacon

%Sustituir valores NAN por ceros.
for i=1:a
    for j=1:INFO{i, 4}
        if isnan(INFO{i, 6}(j, 1)) == true
            INFO{i, 6}(j, 1) = 0;
        end
    end
end 

%Para corregir cuando las señales no terminar o comienzan al mismo tiempo.
% for i=1:81300
%     e_c(i) = INFO{1, 6}(i, 1);
%     %z_c(i) = INFO{3, 6}(i, 1);
% end

% INFO{1, 6} = e_c';
%INFO{3, 6} = z_c';

%Juntar señales cuando hay mas de una hora.
% merge = cell(3, 9);
% if a > 3
%     for i=1:3
%         merge{i,1} =INFO{i,1};
%         merge{i,2} =INFO{i,2};
%         merge{i,3} =INFO{i,3};
%         merge{i,4} =INFO{i,4}+INFO{i+3,4};
%         merge{i,6} =cat(1,INFO{i,6},INFO{i+3,6});
%     end
% end
% INFO=merge;
% a=3;

%Para el caso T4S29 y T4Q99.
% merge = cell(3, 9);
% if a > 3
%     for i=1:3
%         merge{i,1} =INFO{i+(i-1),1};
%         merge{i,2} =INFO{i+(i-1),2};
%         merge{i,3} =INFO{i+(i-1),3};
%         merge{i,4} =INFO{i+(i-1),4}+INFO{i+i,4};
%         merge{i,6} =cat(1,INFO{i+(i-1),6},INFO{i+i,6});
%     end
% end
% INFO=merge;
% a=3;

% for i=1:415600
%     e_c(i) = INFO{1, 6}(i, 1);
%     %z_c(i) = INFO{3, 6}(i, 1);
% end
% 
% INFO{1, 6} = e_c';
% %INFO{3, 6} = z_c';
% 
% INFO{1, 4} = 415600;

%%
%Correcion instrumental. a) Por factores, b) polos y zeros.
corr = 1;

%Multiplicando por factores.
if corr == a

   %Valores ganancia del sensor en [V/m/s]. Orden [EW NS Z]
   gs = [800 800 800];
   %Valores ganancia registrador en [cuentas/V]. Orden [EW NS Z]
   gr = [311915.159 313185.092 314465.409];

   INFO{1, 6} = INFO{1, 6}*(1/gs(1))*(1/gr(1))*1000; %1000 para unidades mm/s
   INFO{2, 6} = INFO{2, 6}*(1/gs(2))*(1/gr(2))*1000;
   INFO{3, 6} = INFO{3, 6}*(1/gs(3))*(1/gr(3))*1000;
else
    %Con polos y zeros
    %pole = [-0.14803+1i*0.14803;-0.14803-1i*0.14803;-314.15927];
    %   -2199.11;-471.24];  %En rad/s
    %pole = [-0.1486+1i*0.1486;-0.1486-1i*0.1486;-1130.97;-1005.31;-502.65];  %En rad/s
    pole = [-0.1486+1i*0.1486;-0.1486-1i*0.1486;-391.95515+1i*850.69303;-391.95515-1i*850.69303;...
        -2199.11486;-471.2389];  %En rad/s
    %pole = [-0.14803+1i*0.14803;-0.14803-1i*0.14803;-314.15927];
    zero = [0.0;0.0];
    %c = [-85435200;-83711400;-84412133.89];
    c = [2.26544e21;2.27776e21;2.25464e21];
    %c = [-7.81572e7;-7.86721e7;-7.89937e7];        %Normalizacion total.

    % n = length(samples);
    % % n = 81300;    %Para señales con tiempos distintos.
    % dt = 1/sps;
    % para.Fmaxdec = 1/dt;
    % f = linspace(0,para.Fmaxdec,n);
    % w = 2*pi*f(1:(n/2+1));

    for i=1:a
        n = length(INFO{i,6});
        dt = 1/INFO{i,2};
        para.Fmaxdec = 1/dt;
        f = linspace(0,para.Fmaxdec,n);
        w = 2*pi*f(1:(n/2+1));
        T = Seismometer_Response(zero,pole,c(i),w);
        sps1 = fft(INFO{i, 6});

        sps1(1:n/2+1) = sps1(1:n/2+1)./T;sps1(1)=0;
        INFO{i, 6} = ifft([sps1(1:n/2+1);conj(sps1(n/2:-1:2))],'symmetric');
        INFO{i, 6} = INFO{i, 6}*1000;    %Para mm/s
        INFO{i, 6} = filtrar(INFO{i,6},1/INFO{i,2},1,100,4,5);
    end
end
%% 
%
%Ploteo de archivos.
%Se realiza el eje de tiempo para los archivos.
for i=1:a
    tiempo=zeros(INFO{i,4},1);
    for j=1:INFO{i,4}
        tiempo(j,1)=INFO{i,3}+(j-1)*(1/(86400*INFO{i,2}));
    end
    INFO{i,5}=tiempo;
end

%Para señales que no estan en tiempo.
% INFO{1, 4} = INFO{2, 4};
% INFO{1, 5} = INFO{2, 5};
% 
% INFO{3, 4} = INFO{1, 4};
% INFO{3, 5} = INFO{1, 5};

for i=1:a
    h=figure(1);
    subplot(a,1,i)
    hold on
    plot(INFO{i,5},INFO{i,6}, 'b', LineWidth=2)
    datetick('x','HH:MM:SS')
    c=INFO{1,3}+1/24;
    %c=INFO{1,3}+2/24;
    xlim([INFO{1,3} c]);
    ax = gca;
    ax.FontSize =13;
    title(['',INFO{i,1}], FontSize=18)
    xlabel('Tiempo, HH:MM:SS')
    ylabel('Velocidad de partícula, mm/s')
    ylim([-10 10])
    grid on
    box on

    xrec = [739153.812313773 739153.856761343 739153.856761343 739153.812313773];
    yrec = [-20 -20 20 20];
    col = [0.8500 0.3250 0.0980];
    fill(xrec, yrec, col, 'FaceAlpha', 0.2)
end
savefig(h,'Señales.fig');

%% 
%Se corta una ventana para analizar el sismo en frecuencias.
[INFOven] = ventana2(INFO); %Se llama a la funcion para cortar el sismo.
save('INFO.mat','INFOven');

%%
%Análisis de Fourier del sismo cortado. Espectro de amplitud.
[INFOFou] = Fourier(INFOven); %Se llama a la funcion para realizar fft.

%Exportar explosion en formato txt.
% carp = uigetdir('C:\Users\erick\Desktop\Proyectos\Procesamiento Huelitli-3D\Explosiones 2da parte\');
% T = table(INFOven{1,6}, INFOven{2, 6}, INFOven{3, 6},...
%     'VariableNames', {'E', 'N', 'Z'});
% writetable(T, fullfile(carp, ['LA11_23' '.txt']), 'Delimiter', ' ');
%% 
%Ploteo de archivos.
tven =50*(INFOven{1, 2}/100);
%Frecuencia vs ppv.
for i = 1:a
    [~, ind] = max(abs(INFOven{i, 6}));
    %[~, ind] = min(abs(info{i, m+2}-info{i, m+5}));
    wind = INFOven{i, 6}(ind-tven:ind+tven);
    fwind = fft(wind)/length(wind);
    esp_wind = abs(fwind);
    esp_wind = esp_wind(1:(length(esp_wind)-1)/2+1);
    esp_wind(2:end-1) = 2*esp_wind(2:end-1);
    INFOven{i, 8} = esp_wind;
    %Eje de frecuencias.
    frec2 = INFOven{i, 2}*(0:(length(wind)/2))/length(wind);
    frec2 = frec2';
    INFOven{i, 7} = frec2;
%     [~, ind2] = max(esp_wind);
%     INFO{i, m+9} = frec(ind2);
end

for i=1:a
    h = figure(4);
    pos1 = [0.05 1-0.3*i 0.6 0.25];
    subplot('Position', pos1)
    plot(INFOven{i, 5}(1:400, 1), INFOven{i, 6}(1:400, 1), 'b', lineWidth=1)
    datetick('x', 'HH:MM:SS')
    title(['', INFOven{i, 1}], FontSize=18)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    ax.FontSize = 13;
    xlabel('Tiempo, HH:MM:SS')
    ylabel('Velocidad de partícula, mm/s')
    %xlim([738789.6825836805 738789.6826824073])
    ylim([-10 10])
    grid on
    [~, indm] = max(abs(INFOven{i, 6}));
    maxt = INFOven{i, 5}(indm);
    maxa = INFOven{i, 6}(indm);
    hold on
    plot(maxt, maxa, '+', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', 'r')
    val1 = INFOven{i, 5}(indm-tven);
    val2 = INFOven{i, 5}(indm+tven);
    xrec = [val1 val2 val2 val1];
    yrec = [-10 -10 10 10];
    col = [0.8500 0.3250 0.0980];
    fill(xrec, yrec, col, 'FaceAlpha', 0.2)
    velmax = round(max(abs(INFOven{i, 6})), 2);
    txt = ['Velocidad de partícula máxima: ' num2str(velmax) ' mm/s'];
    text(INFOven{i, 5}(330),6,txt, 'FontSize', 12, 'FontWeight','bold')
    
    pos2 = [0.7 1-0.3*i 0.25 0.25];
    subplot('Position', pos2)
    plot(INFOven{i, 7}, INFOven{i, 8}, 'b', linewidth=1)
    title(['',INFOFou{i, 1}], FontSize=18)
    ax = gca;
    ax.TitleHorizontalAlignment = 'right';
    ax.FontSize = 13;
    xlabel('Frecuencia, Hz')
    ylabel('Amplitud')
    %set(gca, 'XScale', 'log')
    %set(gca, 'YScale', 'log')
    grid on
    xlim([1 100]);
    ylim([0 0.8])
    [~, indf] = max(abs(INFOven{i, 8}));
    maxf = INFOven{i, 7}(indf);
    maxaf = INFOven{i, 8}(indf);
    hold on
    plot(maxf, maxaf, '+', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', 'r')
    txt = ['Frecuencia máxima: ' num2str(round(maxf)) ' Hz'];
    text(10, 0.7,txt, 'FontSize', 12, 'FontWeight','bold')
end

save('INFO200sps.mat','INFOven');