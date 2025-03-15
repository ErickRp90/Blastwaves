%Procesamiento explosiones.
%Programa para realizar el análisis de vibraciones producidas por
%explosiones.

clearvars
clc

%Lectura de archivos coordenadas de estaciones y pruebas.
%oldfolder = cd('C:\Users\erick\Desktop\Programas\Matlab');
oldfolder = cd('C:\Users\ERamosP\Desktop\');
addpath(oldfolder)
folder = uigetdir('C:\Users\ERamosP\Desktop\\Proyectos\IINGEN\Estudios de vibracion\BLOQUE ALACTE\',...
    'Selecciona carpeta de estaciones y pruebas');
% folder = uigetdir('C:\Users\erick\Desktop\Proyectos\Huelitli-3D\Huelitli 4ta Etapa\',...
%    'Selecciona carpeta de estaciones y pruebas');
files = dir(fullfile(folder, '*.txt'));
nf = length(files);
cd (folder)

fidA = fopen(fullfile(folder, files(1).name));
EstA = textscan(fidA, '%f %s %f %f', 'HeaderLines', 1);
fidB = fopen(fullfile(folder, files(2).name));
EstB = textscan(fidB, '%f %s %f %f', 'HeaderLines', 1);
fidPA = fopen(fullfile(folder, files(3).name));
PruA = textscan(fidPA, '%f %f %f %f %f', 'HeaderLines', 1);
fidPB = fopen(fullfile(folder, files(4).name));
PruB = textscan(fidPB, '%f %f %f %f %f', 'HeaderLines', 1);

%Cargar archivos localización A.
loca_fold = uigetdir('C:\Users\ERamosP\Desktop\Proyectos\IINGEN\Estudios de vibracion\BLOQUE ALACTE\',...
    'Selecciona carpeta localización');
% loca_fold = uigetdir('C:\Users\erick\Desktop\Proyectos\Huelitli-3D\Huelitli 4ta Etapa\',...
%    'Selecciona carpeta localización B');
loca_files = dir(fullfile(loca_fold, '*.txt'));
nafiles = length(loca_files);

nexp = 3;%24;
nesta = nafiles/nexp;
info = cell(nafiles, 12);
%Procesamiento para grafica de distancia vs ppv y
%grafica frecuencia vs ppv
for i = 1:nafiles
    codea = split(loca_files(i).name, [".", "_"]);
    info{i, 1} = string(codea(1));  %Nombre
    info{i, 2} = str2double(codea(2));  %# explosion
    fid = fopen(fullfile(loca_fold, loca_files(i).name));
    expdata = textscan(fid, '%f %f %f', 'HeaderLines', 1);
    expdata = cell2mat(expdata);
    info{i, 3} = filtrar(expdata(:,1), 1/100, 1, 50, 4, 5);    %EW
    info{i, 4} = filtrar(expdata(:,2), 1/100, 1, 50, 4, 5);    %NS
    info{i, 5} = filtrar(expdata(:,3), 1/100, 1, 50, 4, 5);    %Z

    %if
    %estm = "LA05";
    %estm = "PA03";
%     if strcmp(info{i, 1}, estm)
%         info{i, 3} = 3 * info{i, 3};
%         info{i, 4} = 3 * info{i, 4};
%         info{i, 5} = 3 * info{i, 5};
%     end

    info{i, 6} = max(abs(info{i,3}));    %ppv EW
    info{i, 7} = max(abs(info{i,4}));    %ppv NS
    info{i, 8} = max(abs(info{i,5}));    %ppv Z
    
    
    %Calculo de la distancia.
    for j = 1:nexp
        if info{i, 2} == PruA{1, 1}(j, 1)
            xp = PruA{1, 2}(j, 1);
            yp = PruA{1, 3}(j, 1);
        end
    end
    for k = 1:nesta
        if strcmp(info{i, 1}, EstA{1, 2}(k, 1))
            xe = EstA{1, 3}(k, 1);
            ye = EstA{1, 4}(k, 1);
        end
    end

    info{i, 9} = sqrt((xp-xe)^2 + (yp-ye)^2);

    %Frecuencia vs ppv.
    for m = 1:3
        [~, ind] = max(info{i, m+2});
        %[~, ind] = min(abs(info{i, m+2}-info{i, m+5}));
        wind = info{i, m+2}(ind-50:ind+50);
        fwind = fft(wind)/length(wind);
        esp_wind = abs(fwind);
        esp_wind = esp_wind(1:(length(esp_wind)-1)/2+1);
        esp_wind(2:end-1) = 2*esp_wind(2:end-1);

        %Eje de frecuencias.
        frec = 100*(0:(length(wind)/2))/length(wind);
        frec = frec';

        [~, ind2] = max(esp_wind);
        info{i, m+9} = frec(ind2);
    end
end


%Graficas.
%%
%Graficas.

%Tiempo.
dt = 0.01;
for i = 1:nafiles
    t = (-1:dt:dt*(length(info{i, 3})-1)-1)';
    info{i, 13} = t;
end

%Grafica norma NOM-026-SESH-2007
norma(:, 1) = (10:10:120);
norma(:, 2) = [19.5, 13, 9, 7, 5.5, 5, 4, 3, 2, 1.5, 1, 0.8];

fnorm = [1;4;16;40;100];
vnorm = [2.54;9;9;25;25];
cmp1 = hsv(nesta);
cmp2 = zeros(3*length(cmp1), 3);
for i=1:length(cmp1)
    cmp2(i+2*(i-1), 1:3) = cmp1(i, 1:3);
    cmp2(i+(2*i-1), 1:3) = cmp1(i, 1:3);
    cmp2(3*i, 1:3) = cmp1(i, 1:3);
end
set(0,'DefaultLegendAutoUpdate','off')
for i = 1:nexp
    h = figure(1);
    sgtitle(sprintf('Localización 1A. Prueba %d. Profundidad = %d m. Carga = %d kg',...
             PruA{1, 1}(i), PruA{1, 4}(i,1), PruA{1, 5}(i,1)), 'FontSize', 20, 'FontWeight', 'bold')
    hold on
    for j = 1:nesta
        axpe = subplot("Position", [0.05, 0.1, 0.20, 0.80]);
        hold on
        colororder(axpe, cmp1);
        plot(info{i+nexp*(j-1), 13}, info{i+nexp*(j-1), 3}+10*(nesta-j), 'LineWidth', 2)
        title('Componente Transversal', 'FontSize', 16, 'FontWeight', 'bold')
        xlabel('Tiempo, s', 'FontSize', 16)
        ylabel('Velocidad de partícula, mm/s', 'FontSize', 16)
        xlim([-1 4])
        ylim([-5 10*nesta])
        grid on
        box on
        ax= gca;
        %ax.FontSize = 12;
        ax.LineWidth = 1;
        ax.YTick = -5:5:10*nesta;
        
        legend_labels{j} = sprintf('%s - %3.0f [m]', info{i+nexp*(j-1), 1}, info{i+nexp*(j-1), 9});

        [~, indm] = max(abs(info{i+nexp*(j-1), 3}));
        maxe(j, 1) = t(indm);
        maxe(j, 2) = info{i+nexp*(j-1), 3}(indm)+10*(nesta-j);
    end
    plot(axpe, maxe(:, 1), maxe(:, 2), '+', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', 'r')
    legend(legend_labels{:}, 'position', [0.75 0.05 0.15 0.15], 'FontSize', 12, 'FontWeight', 'bold', 'NumColumns', 2);

    for j = 1:nesta
        axpn = subplot("Position", [0.25, 0.1, 0.20, 0.80]);
        hold on
        colororder(axpn, cmp1);
        plot(info{i+nexp*(j-1), 13}, info{i+nexp*(j-1), 4}+10*(nesta-j), 'LineWidth', 2)
        title('Componente Radial', 'FontSize', 16, 'FontWeight', 'bold')
        xlabel('Tiempo, s', 'FontSize', 16)
        xlim([-1 4])
        ylim([-5 10*nesta])
        yticks([])
        grid on
        box on
        ax= gca;
        %ax.FontSize = 12;
        ax.LineWidth = 1;
        ax.YTick = -5:5:10*nesta;
        ax.YTickLabel = [];

        [~, indm] = max(abs(info{i+nexp*(j-1), 4}));
        maxe(j, 1) = t(indm);
        maxe(j, 2) = info{i+nexp*(j-1), 4}(indm)+10*(nesta-j);
    end
    plot(axpn, maxe(:, 1), maxe(:, 2), '+', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', 'r')

    for j = 1:nesta
        axpn = subplot("Position", [0.45, 0.1, 0.20, 0.80]);
        hold on
        colororder(axpn, cmp1);
        plot(info{i+nexp*(j-1), 13}, info{i+nexp*(j-1), 5}+10*(nesta-j), 'LineWidth', 2)
        title('Componente Vertical', 'FontSize', 16, 'FontWeight', 'bold')
        xlabel('Tiempo, s', 'FontSize', 16)
        xlim([-1 4])
        ylim([-5 10*nesta])
        yticks([])
        grid on
        box on
        ax= gca;
        %ax.FontSize = 12;
        ax.LineWidth = 1;
        ax.YTick = -5:5:10*nesta;
        ax.YAxisLocation = 'right';

        [~, indm] = max(abs(info{i+nexp*(j-1), 5}));
        maxe(j, 1) = t(indm);
        maxe(j, 2) = info{i+nexp*(j-1), 5}(indm)+10*(nesta-j);
    end
    plot(axpn, maxe(:, 1), maxe(:, 2), '+', 'MarkerSize', 12, 'LineWidth', 1.5, 'Color', 'r')        
    
    for j = 1:nesta
        axdis = subplot("Position", [0.7, 0.60, 0.25, 0.30]);
        hold on
        colororder(axdis, cmp2);
        scatter(info{i+nexp*(j-1), 9}, info{i+nexp*(j-1), 6}, 70, 'filled')
        scatter(info{i+nexp*(j-1), 9}, info{i+nexp*(j-1), 7}, 70, 'square', 'filled')
        scatter(info{i+nexp*(j-1), 9}, info{i+nexp*(j-1), 8}, 70, 'diamond', 'filled')
        if j == nesta
            plot(norma(:, 1), norma(:, 2), 'Color', 'k', LineWidth=2)
            xline(300,'--r',{'300 m'}, 'Linewidth', 2, 'Fontsize', 14);
        end
        xlabel('Distancia, m', 'FontSize', 14)
        ylabel('Velocidad de partícula, mm/s', 'FontSize', 14)
        xlim([1 500])
        ylim([0.1 100])
        %set(gca, 'Xscale', 'log')
        set(gca, 'Yscale', 'log')
        grid on
        box on
        ax= gca;
        %ax.FontSize = 12;
        ax.LineWidth = 1;
        %ax.XTickLabel = [1 100 500];
        ax.YTickLabel = [0.1 1 10 100];
        if j == 1
            legend('Transversal', 'Radial', 'Vertical', 'Location','southwest', 'FontSize', 15, 'Fontweight', 'bold')
        end

        axfre = subplot("Position", [0.7, 0.25, 0.25, 0.30]);
        hold on
        colororder(axfre, cmp2);
        scatter(info{i+nexp*(j-1), 10}, info{i+nexp*(j-1), 6}, 70, 'filled')
        scatter(info{i+nexp*(j-1), 11}, info{i+nexp*(j-1), 7}, 70, 'square', 'filled')
        scatter(info{i+nexp*(j-1), 12}, info{i+nexp*(j-1), 8}, 70, 'diamond', 'filled')
        if j == nesta
            plot(fnorm, vnorm, 'Color', 'k', LineWidth=2)
        end
        xlabel('Frecuencia, Hz', 'FontSize', 14)
        ylabel('Velocidad de partícula, mm/s', 'FontSize', 14)
        xlim([1 50])
        ylim([0.1 100])
        set(gca, 'Xscale', 'log')
        set(gca, 'Yscale', 'log')
        grid on
        box on
        ax= gca;
        %ax.FontSize = 12;
        ax.LineWidth = 1;
        ax.XTickLabel = [1 10 50];
        ax.YTickLabel = [0.1 1 10 100];
    end
    set(h, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'white')
    saveas(h, sprintf('%d - Prueba 1A', i), 'jpg')
    saveas(h, sprintf('%d - Prueba 1A', i), 'fig')
    close(h)
end

fclose('all');