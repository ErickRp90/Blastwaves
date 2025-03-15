function xc = corre_ins(x, sps, p, z, c, tipo)
%Funcion para realizar la correccion instrumental con Polos y zeros.
%       x: senal a corregir.
%       sps: samples per second.
%       p: vector con polos, Hz.
%       z: vector con zeros, Hz.
%       c: total normalization factor.
%
%Tipo:  1) Displacement.
%       2) Velocity.
%       3) Acceleration.

%%
%Transformada de Fourier de la senal. Se trabaja hasta la frecuencia de Nyquist.
samples = length(x);
% samples = 2^nextpow2(samples);          %Mejorar desempeno de fft.
fftx = fft(x, samples);
fasx = abs(fftx/samples);
fasx1 = fasx(1:samples/2);          %Se considera muestras pares.
freq = (sps/samples)*(0:samples-1)'; %Entre 2
ome = 2*pi*freq;

%%
%Funcion de trasferencia con los polos y zeros.
%T = C* Z/P
np = length(p);
nz = length(z);
if tipo == 1
    num = (1i*ome-z(1)).*(1i*ome-z(2)).*(1i*ome-z(3)); %No hay z(3)
elseif tipo == 2
    num = (1i*ome-z(1)).*(1i*ome-z(2));
elseif tipo == 3
    num = (1i*ome-z(1));
end
den = (1i*ome-p(1));
for i =2:np
    den = den.*(1i*ome-p(i));
end

TF = c*(num./den);
TF = ((1i*ome-z(1)).*(1i*ome-z(2)))./((1i*ome-p(1)).*(1i*ome-p(2)).*(1i*ome-p(3)).*(1i*ome-p(4))...
    .*(1i*ome-p(5)).*(1i*ome-p(6)));
TF = c*TF;
length(TF)
TF(1,1) = TF(2,1);        %Ya que los zeros son iguales a 0.
%TF_fas = abs(TF);
TF_fas = sqrt(real(TF).^2+imag(TF).^2);
%Pha = phase(TF);
Pha = atand(imag(TF)./real(TF));

%%
%Graficas.
figure(1)
hold on
%plot(fre, ASx1, 'b')
plot(ome/(2*pi), TF_fas, 'r')
%plot(fre, 1./Amp, 'k')
plot(1,1,'o')
grid on
xlim([0.01 200])
xlabel('Frecuencia, Hz')
ylabel('Amplitud, cuentas/m ')
ylim([1e-2 1e1])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

figure(2)
plot(ome/(2*pi), 20*log10(TF_fas), 'r')
grid on
xlim([0.01 200])
xlabel('Frecuencia, Hz')
ylabel('Amplitud, dB ')
ylim([-40 20])
set(gca, 'XScale', 'log')

figure(2)
hold on
%plot(fre, Pha*(180/pi), 'b')
grid on
% xlim([0.01 50])
ylim([-100 100])
xlabel('Frecuencia, Hz')
ylabel('Fase, grados')
set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')

%%
%Senal corregida.
XT = FT_x(1:samples/2+1)./TF;
size(XT)
size(TF)
XTco = conj(XT(2:end-1));
XTco2 = fliplr(XTco);
XTC = cat(1,XT,XTco2);
size(XTC)
xc = ifft(XT, 'symmetric');
length(xc)

%Filtrar la señal 
fs = 100;
fn = fs/2;
lowcut = 0.1/fn;
highcut = 50/fn;
order = 4;

%Design the Butterworth filter.
[b, a] = butter(order, [lowcut, highcut]/(fs/2), 'bandpass');
xc = filtfilt(b, a, xc);

XCF = fft(xc);
AXCF = abs(XCF/samples);
AXCF2 = AXCF(1:samples/2+1);
plot(fre, AXCF2, 'r')
% xc = real(xc);
% t = (0:1:samples-1)';
figure(3)
plot(xc)