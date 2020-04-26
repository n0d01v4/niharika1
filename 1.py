# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 15:00:11 2019

@author: NIHARIKAPRIYADARSHI
"""
from array import *
import numpy as np
import math
#Input parameters
#Cambered vane for specified diameter ratio
R=287
G=1.323
G1=0.5*(G-1)
G2=G/(G-1)
G3=1/G2
G4=-(G+1)/(2*(G-1))
pi=3.14159
Rad=pi/180
alpha0=55.6
alpha1=72
T01=1200       #Turbine inlet stagnaion temperature
P=64.5*1000    #Power in watt
w=0.22         #mass flow rate in kg/s
Cp=1184
N=118000       #speed in rpm
P00=385*1000   #kpa
RHO00=P00/(R*T01)
R1R1A=1.25     #stator exit to rotor inlet pressure ratio
RT2R1A=0.7133  #rotor exit tip to rotor inlet radius ratio
R0R1A=1.2581   #Stator inlet to rotor inlet radius ratio
RH2RT2=0.3493  #Rotor exit hub to tip radius ratio
CDT2=0.0023    #ratio of clearance gap
MU=0.5804*0.0001        #Gas viscosity
K=5            #number of radial sectors at rotor exit
LCDH=CDT2*2/(1-RH2RT2)
LDF=0.0         #L_df
DHSHFT=P/w     #Delta shaft
DHVDAV=(DHSHFT+LDF)/(1-LCDH)

Nr=int((pi/30)*(110-alpha1)*math.tan(alpha1*Rad))     # Number of blades
x=1-(2/Nr)
RV1AAV=1
while DHVDAV>=0:
    U1A=(DHVDAV*RV1AAV/x)**0.5
    D1A=(U1A*2*60)/(2*pi*N)
    VU1A=U1A*x
    
    #Stator exit
    D1=R1R1A*D1A
    VU1=VU1A/R1R1A
    es=0.05
    V1=VU1/math.sin(alpha1*Rad)
    VR1=VU1/math.tan(alpha1*Rad)
    T1=T01-((V1**2)/(2*Cp))
    V1IDS=(V1**2)/(1-es)
    P1=(1-((V1IDS)/(2*Cp*T01)))*P00
    P01=P1*(T01/T1)**G2
    RHO1=P1/(R*T1)
    Hs=w/(RHO1*VR1*pi*D1)  #Stator height
    
    Vcr1=(2*G*R*T01/(G+1))**0.5
    
    #Stator inlet
    cs=0.5          #Stator solidity
    D0=D1A*R0R1A
    FICM=0.5*(alpha0+alpha1)                  #Stagger angle
    q=((D0**2-D1**2)/math.cos(FICM*Rad))**2
    r=(((D0**2+D1**2)**2)-q)**0.5
    Chord=(((D0**2+D1**2)-r)**0.5)/2
    hcs=Hs/Chord         #(h/c)_s
    #0cam=alpha1-alpha0-((math.acos(((D0**2)+(D1**2)-(2*Chord)**2)/(2*D0*D1)))/Rad)       #Camber angle
    ss=Chord/cs
    ns=int(pi*D1/ss)
    ##Needs feedback mechanism to be implemented
    while RHO01<=RHO00:
        V0=w/(pi*D0*Hs*RHO01*math.cos(alpha0*Rad))
        T0=T01-((V0**2)/(2*Cp))
        P0=P00*(T0/T01)**G2
        RHO0=P0/(R*T0)
        if RHO0==RHO01:
            break
        RHO01=RHO0
    VU0=V0*math.sin(alpha0*Rad)
    VR0=V0*math.cos(alpha0*Rad)       #Radial velocity
    V0VCR0=V0/Vcr1
    #lcs=pi*0cam/(360*math.sin(0cam*0.5*Rad))
    
    #Rotor inlet
    V1A=VU1A/math.sin(alpha1*Rad)
    T1A=T01-((V1A**2)/(2*Cp))
    P1A=P01*(T1A/T01)**G2
    RHO1A=P1A/(R*T1A)
    Vrad1A=VR1*RHO1*D1/(RHO1A*D1A)
    alpha1A=(math.atan(VU1A/Vrad1A))/Rad
    LDFf=(0.02125*RHO1A*(U1A**3)*(0.5*D1A)**2)/(w*(RHO1A*U1A*0.5*D1A/MU)**0.2)
    DHVDAVf=(DHSHFT+LDFf)/(1-LCDH)
    if DHVDAVf==DHVDAV:
        DHVDAV=DHVDAVf
        break
    DHVDAV=DHVDAVf
WU1A=VU1A-U1A
Beta1A=(math.atan(WU1A/Vrad1A))/Rad
W1A=Vrad1A/math.cos(Beta1A*Rad)
T01AR=T1A+((W1A**2)/(2*Cp))
P01AR=P1A*(T01AR/T1A)**G2
Wcr1A=Vcr1*(T01AR/T01)**0.5
VOVCRA=V1A/Vcr1
WOWCRA=W1A/Wcr1A
print(Beta1A)

# rotor exit
VU2M=200           #Need to specify  vu2m
KK=K+2
D2(KK)=D1A*RT2RIA  #Rotor exit tip diameter
i=1
while i<=kk:
    VX2M=230
    D2.append(D2(KK)*i/(k+2))
    U2.append(pi*D2[i]*N/(2*pi*60))
    T2IR.append(T01AR+((U2[i]**2-U1A**2)/2*Cp))
    Wcr2I.append(Wcr1A*(T2IR/T01AR))
    P2IDI.apppend(P01AR*(T2IR/T01AR)**G2)
    i=i+1
er=             #rotor loss coefficient 
WU2M=VU2M-U2[4]
Beta2M=(math.atan(WU2M/VX2M))/Rad         #VX2M NOT DECLARED YET
W2M=VX2M/math.cos(Beta2M*Rad)
T2M=T2IR[4]-((W2M**2)/(2*Cp))
W2IDM=(W2M**2)/(1-er)
P2M=((1-((W2IDM**2)/(2*Cp*T2MR)))**G2)*P2IDM
RHO2M=P2M/(R*T2M)

while i<=kk:
    
    VU2J.append(VU2M*RV2j2m*D2M/D2[i])   ###RV2j2m needs to be specified
    WU2J.append(VU2J-U2[i])
    Pi=P2M
    RHOi=RHO2M
    if i>(int(kk*0.5)+1):
        P2J.append(Pi+RHOi*0.5*(((VU2J[i]**2)/D2[i])+((VU2J[i+1]**2)/D2[i+1]))*(D2[i]-D2[i+1]))
        W2IDJ.append((2*Cp*T2IR[i]*(1-(P2J[i]/P2IDI[i])**(1/G2)))**0.5)
        W2J.append(((W2IDJ[i]**2)*(1-er))**0.5)
        Beta2j.append((math.asin(WU2J[i]/W2J[i]))/Rad)
        VX2J.append(W2J[i]*math.cos(Beta2j[i]/Rad))
        T2J.append(T2IR[i]-((W2J[i]**2)/(2*Cp)))
        RHO2J.append(P2J[i]/(R*T2J[i]))
        Pi=P2J[i]
        RHOi=RHO2J[i]
    elif i<(int(kk*0.5)+1):
        P2J.append(Pi+RHOi*0.5*(((VU2J[i]**2)/D2[i])+((VU2J[i-1]**2)/D2[i-1]))*(D2[i]-D2[i-1]))
        W2IDJ.append((2*Cp*T2IR[i]*(1-(P2J[i]/P2IDI[i])**(1/G2)))**0.5)
        W2J.append(((W2IDJ[i]**2)*(1-er))**0.5)
        Beta2j.append((math.asin(WU2J[i]/W2J[i]))/Rad)
        VX2J.append(W2J[i]*math.cos(Beta2j[i]/Rad))
        T2J.append(T2IR[i]-((W2J[i]**2)/(2*Cp)))
        RHO2J.append(P2J[i]/(R*T2J[i]))
        Pi=P2J[i]
        RHOi=RHO2J[i]
    elif i==(int(kk*0.5)+1):
        P2J.insert(int(kk*0.5)+1,P2M)
        W2IDJ.insert(int(kk*0.5)+1,W2IDM)
        W2J.insert(int(kk*0.5)+1,W2M)
        Beta2j.insert(int(kk*0.5)+1,Beta2M)
        VX2J.insert(int(kk*0.5)+1,VX2M)
        T2J.insert(int(kk*0.5)+1,T2M)
        RHO2J.insert(int(kk*0.5)+1,RHO2M)
    i=i+1
        
t=0
while i<=kk:
    t=t+RHO2J[i]*VX2J[i]*D2[i]
    i=i+1
w_cal=t*pi*(D2t-D2h)/(2*K)

while i<=kk
DHVDI.append((U1A*VU1A)-(U2[i]*VU2J[i]))
T02I.append(T0-(DHVDI[i]/Cp))
Vcr2I.append(Vcr1A*(T02I[i]/T01)**0.5)
alpha2I.append((math.atan(VU2J[i]/VX2J[i]))/Rad)
P02I.append(P2J[i]*(T02I[i]/T2J[i])**G2)
V2I.append(VX2I/math.cos(alpha2I[i]*Rad))
i=i+1

#performance
Ls=(V1id**2-V1**2)    ##need o search what v1iid is
Lrav=0
Lexav=0
while i<=kk:
    Lrav=Lrav+((wi/w)*(W2idi**2-W2i**2)/2)
    Lexav=Lexav+((wi/w)*(V2i**2)/2)
    i=i+1
    DH0IDI.append(Cp*T01*(1-(P02I[i]/P01)**(1/G2)))
    DHIDI.append(Cp*T01*(1-(P2J[i]/P01)**(1/G2)))
    n0VDi=DHVDI[i]/DH0IDI[i]
    nVDi=DHVDI[i]/DHIDI[i]
    i=i+1
#overall veocity diagram efficiency
t=0
n=0
while i<=kk:
    t=t+(wi/w)*DH0IDI
    n=n+(wi/w)*DHIDI
    i=i+1
n0VDAV=DHVDAV/t
nVDAV=DHVDAV/n
#net overall efficiency
n0shft=n0VDAV*DHshaft/DHVDAV
nshft=nVDAV*DHshaft/DHVDAV
DHIDAV=0
while i<=kk:
    DHIDAV=DHIDAV+DH0IDI[i]
    i=i+1
DHIDAV=DHIDAV/KK
N_sp=(N*(w/RHO2)**0.5)/(DHIDAV)**(3/4)           #Specific speed
