function [ftran,ome,Nome]=FFTS(F,omega,fs,s,Sigtype)

NF=length(F); NFR=round(NF/2);
Freqs_y=fft(F,NF);
if strcmp(Sigtype,'G')
    Pyys_y=(Freqs_y)/fs;
elseif strcmp(Sigtype,'P')
    Pyys_y=(Freqs_y)/NFR;%Freq1.*conj(Freq1)/N;
end

if s==0
    ome=omega; 
    Ns=NF/fs;
    Nome=round(ome*Ns/(2*pi))+1;
    if ome==0
        ftran=Pyys_y(Nome)/2;
    else
        ftran=Pyys_y(Nome);
    end
elseif s==1
    ome=2*pi*(0:NFR)*fs/NF;%f=rad/s
    Nome=[1:NFR+1];
    Pyys_y(1)=Pyys_y(1)/2;
    ftran=Pyys_y(1:NFR+1);
else
   error('F(iw)-0;F_scan data-1;');
end
end
