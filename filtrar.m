function acf = filtrar(ac,dt,Fmin,Fmax,Npol,tap)
%FILTRAR  Aplica tapering a una señal, calcula EAF, filtra las frecuencias,
%  aplica la antitransformada de Fourier y recupera la señal filtrada.
No=length(ac);
Np=2^(fix(log(No)/log(2))+1);% inicio del tapering.
Ntap=fix(tap/100*No);
f_tap(1:Ntap)=sin((0:Ntap-1)*(pi/2)/Ntap);
f_tap=f_tap';
ac(1:Ntap)=ac(1:Ntap).*f_tap(1:Ntap);
ac(No:-1:No-Ntap+1)=ac(No:-1:No-Ntap+1).*f_tap(1:Ntap);
ac(No+1:Np)=0;
Awc=fft(ac,Np);
%sum_ac=Awc(1);
Awc(1)=[];
Awc(Np/2+1:Np-1)=[];
Fnq=1/(2*dt);
%Df=1/(2*Np/2*dt);
Fr=((1:Np/2)/(Np/2)*Fnq)';
%Fil=1./((1+(Fr./Fmin).^(2*Npol)).^0.5);%pasabaja
%Fil=1./((1+(Fmin/Fr).^(2*Npol)).^0.5);%pasalta
Fil=1./sqrt(1.+((Fr.^2.-Fmin*Fmax)./(Fr.*(Fmax-Fmin))).^(2*Npol)); %pasabanda
%Fil=1./((1+(Fr./Fmin).^(2*Npol)).^0.5)+1./((1+(Fmax./Fr).^(2*Npol)).^0.5); %detenbanda
Awcf=Awc.*Fil;
%Awc2(1,1)=sum_ac;
Awc2(1,1)=0;
Awc2(2:Np/2+1,1)=Awcf(1:Np/2);
Awc2(Np/2+2:Np,1)=conj(Awcf(Np/2-1:-1:1,1));
acf=real(ifft(Awc2,Np));
if Np>No
acf(No+1:Np)=[];
end
acf=acf(1:No);