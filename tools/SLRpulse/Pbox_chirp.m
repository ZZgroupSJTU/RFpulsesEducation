function Pbox_chirp(pbw,pd,offset,dt,flip,tpwr,pw,B1map,adb)
%%
% pbw:   pulse band width (Hz)
% pd:    pulse duration   (s)
% offset:Pulse offset  (Hz)
% dt:    The dewell tme 
% flip1: The beginning flip angle of the chirp pulse (TE weighting correction or enrst angle adjust) 
% flip2: The end flip angle of the chirp pulse (TE weighting correction or enrst angle adjust) 
%        When flip2=flip1 means using the same flip angle
% tpwr:  Reference power of hard pulse
% pw:    Reference duration of hard pulse
% adb:   factor of power 
%%
chirp.pbw=20000;
chirp.pd=0.01;
chirp.offset=0;
chirp.dt=4e-6;
chirp.flip=90;
chirp.reftpwr=54;
chirp.refpw=20e-6;
chirp.B1map=1;
chirp.adb=4;

if nargin >1
    chirp.pbw=pbw;
    if nargin>2
        chirp.pd=pd;
        if nargin>3
            chirp.offset=offset;
            if nargin>4
                chirp.dt=dt;
                if nargin>5;
                    chirp.flip=flip;
                    if nargin>6
                        chirp.reftpwr=tpwr;
                        if nargin>7
                            chirp.refpw=pw;
                            if nargin>8
                                chirp.B1map=B1map;
                                if nargin>9
                                    chirp.adb=adb;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if chirp.flip<=pi/2
    FactorB1max=0.26*chirp.flip/(pi/2);
elseif abs(chirp.flip-pi)<pi/6
    FactorB1max=0.26*3;
else
    disp('pay attention - not small flip angle');
end

Oi=-chirp.pbw/2-chirp.offset;
Of=chirp.pbw/2-chirp.offset;
R=(Of-Oi)/chirp.pd;

B1max=1e-3*FactorB1max*sqrt(abs(R));
B1maxref=0.25*1e-3/pw;
rfpwr= refpwr+20*log10(B1max/B1maxref);
rfpwr1=4095*10^(-(ceil(rfpwr)-rfpwr)/20);
chirp.PowerCoarse=ceil(rfpwr);
chirp.PowerFine=rfpwr1;

chirp.t=0:chirp.dt:chirp.pd-chirp.dt;
% if length(chirp.B1map)>1
%     chirp.B1map=interp1(linspace(-2,2,length(chirp.B1map)),1./chirp.B1map,0.+linspace(-2,2,length(t)),'linear','extrap');

chirp.B1amp=(1-(cos(pi*chirp.t/chirp.pd)).^80).*chirp.B1map;
chirp.B1phase=2*pi*(Oi*t+0.5*R*(chirp.t).^2);
B1amp=1023*chirp.B1/max(chirp.B1);

B1phase=

