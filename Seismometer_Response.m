function T=Seismometer_Response(ZZero,pole,gain,w)

if ~iscolumn(ZZero)
    ZZero=ZZero.';
end
if ~iscolumn(pole)
    pole=pole.';
end

n=length(w);
T=zeros(n,1);
for i=1:n
    T(i)   = gain*prod(w(i)-ZZero)/prod(w(i)-pole);
    %T(i)   = gain*prod(1i*w(i)-ZZero)/prod(1i*w(i)-pole);
end
% Amp = gain*abs(T);
% Phi = atan(imag(T)./real(T));
% 
% figure;
% subplot(2,1,1)
% plot(w/2/pi, 20*log10(Amp/Amp(1000)),'linewidth',3)
% ylabel('Amplitud','fontsize',16')
% xlabel('Frecuencia (Hz)','fontsize',16')
% set(gca,'fontsize',16','xscale','log')
% set(gca, 'YScale', 'log')
% xlim([0.01 100])
% 
% subplot(2,1,2)
% plot(w/2/pi,Phi,'linewidth',3)
% ylabel('Fase','fontsize',16')
% xlabel('Frecuencia (Hz)','fontsize',16')
% set(gca,'ytick',[-pi/2;0;pi/2],'yticklabel',['-pi/2';'    0';' pi/2'],'fontsize',16','xscale','log')
% xlim([0.01 100])
% set(gcf,'color',[1 1 1])
