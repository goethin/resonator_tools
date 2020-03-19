from resonator_tools import circuit
import numpy as np
import matplotlib.pyplot as plt

port1 = circuit.notch_port_lhalfhangers()

#port1.add_froms2p('S21testdata.s2p',3,4,'realimag',fdata_unit=1e9,delimiter=None)
port1.add_froms2p('hangers.s2p',1,2,'dBmagphasedeg')
#port1.add_froms2p('C:\Users\Jann\switchdrive\Resonators\Data\Bluefors\NbTiN20\NbTiN20_Lhalfhangers\hangerm30dB100Hz22mK.csv',1,2,'dBmagphaserad',1,delimiter=None)
port1.cut_data(3.975e9,4e9)
port1.autofit()
print("Fit results:", port1.fitresults)
port1.plotall()
print("single photon limit:", port1.get_single_photon_limit(diacorr=True), "dBm")
print("photons in reso for input -140dBm:", port1.get_photons_in_resonator(-150,unit='dBm',diacorr=True), "photons")
print("done")

results=port1.fitresults


def convert_fitresults_to_lambda_half_hanger(fitresults):
    """
    PREVIOUS CODE:
    The code of [0], fits a circle in the (S21_real,S21_imag) space to the resonance of a quareter-wavelength notch-type resonator.
    He distinguishes between the approximative PhiRM method and the diameter corrected method.    [1] introduces the PhiRM method and [2] differentiates from the correct method.
    The program fits the circle with the algebraic technique described in [3]

    [0] (Main paper) Probst, et al."Efficient and robust analysis of complex scattering data under noise in microwave resonators" Review of Scientific Instruments, 86(2), p.024706 (2015).
    [1] (not diameter corrected) Jiansong Gao, "The Physics of Superconducting Microwave Resonators" (PhD Thesis), Appendix E, California Institute of Technology, (2008)
    [2] (diameter corrected) M. S. Khalil, et. al., J. Appl. Phys. 111, 054510 (2012)
    [3] (fitting techniques) N. CHERNOV AND C. LESORT, "Least Squares Fitting of Circles", Journal of Mathematical Imaging and Vision 23, 239, (2005)
    [4] (further fitting techniques) P. J. Petersan, S. M. Anlage, J. Appl. Phys, 84, 3392 (1998)
    """
    """
    NEW (implemented by Jann Ungerer, March 17, 2020):
    If the notch-type resonator (hanger) is a half-wavelength resonator, the definition of Qc changes by a factor of 2
    (See e.g. activity report of Jann from November 2019.)
    When using the correct fitting method (diameter corrected, described by [2]),
    as a result also the total Ql and Qi have to be adapted.
    We use the framework of [0] for 
    """
    r=fitresults
    Qe_lhalf=2*r['Qc_dia_corr']
    Qi_lhalf=1/(1/r['Ql']-1/r['Qc_dia_corr'])
    Ql_lhalf1=1/(1/Qi_lhalf+1/Qe_lhalf)
   
    print("lhalfhangers:\nQi=",Qi_lhalf,"\nQe=",Qe_lhalf,"\nQl=",Ql_lhalf1,"\nQl2=",Ql_lhalf2)
    pass
    

QLold=results['Ql']
Qiold=results['Qi_no_corr']
Qcold=results['absQc']
phi0=results['phi0']
fr=results['fr']


def S21_notch(f,fr,Ql,Qc,phi,a=1.,alpha=0.,delay=0.):
    '''
    full model for notch type resonances
    '''
    return a*np.exp(np.complex(0,alpha))*np.exp(-2j*np.pi*f*delay)*(1.-Ql/Qc*np.exp(1j*phi)/(1.+2j*Ql*(f-fr)/fr)) 


def S21_notch_lhalf(f,a,f0,Qe,Qi,phi,delay):
    Qe=Qe/2#Difference in fitting formula for half-wavelength resonator compared to quarter-wavelength resonator. It comes from the definition of Qc as in Thomas Hasler's phd thesis (p.115&116)
    Ql=1/(1/Qe+1/Qi)
    alpha=0
    Qe=np.abs(Qe)#Q>0
    Qi=np.abs(Qi)
    S21=a*np.exp(1j*alpha)*np.exp(-2*np.pi*1j*f*delay)*(1-(np.exp(1j*phi)*Ql/Qe)/(1+2j*Ql*(f/f0-1)))
    return np.abs(S21)

f=port1.f_data


def absolute(data):
    real=getattr(data,"real")
    imag=getattr(data,"imag")
    return np.sqrt(real**2+imag**2)

r=results
radius=r['Ql']/r['Qc_dia_corr']
Qcnew=2*r['Qc_dia_corr']
Qlnew=2956.5535186146076 
Qinew=1/(1/Qlnew-1/Qcnew)
Qinew2=r['Qi_dia_corr']
Qlnew2=1/(1/Qinew2+1/Qcnew)
kappanew=fr/Qlnew
kappanew2=fr/Qlnew2

kappaold=fr/r['Ql']

print("Q factors for half-wavelength resonator",Qlnew,Qcnew,Qinew)

plt.figure()
plt.plot(f,absolute(port1.z_data),label="zdata")
plt.plot(f,np.abs(S21_notch(f,fr,Ql=QLold,Qc=Qcold,phi=phi0,a=1.,alpha=0.,delay=0.)),label="notch fit")

plt.plot(f,np.abs(S21_notch_lhalf(f,a=1,f0=fr,Qe=2*Qcold,Qi=Qiold,phi=phi0,delay=0.)),label="l half fit")
plt.legend()

f=np.linspace(3.98e9,3.99e9,1e6)
S21min=min(np.abs(S21_notch(f,fr,Ql=QLold,Qc=Qcold,phi=phi0,a=1.,alpha=0.,delay=0.)))

fmin=f[np.argwhere(S21_notch(f,fr,Ql=QLold,Qc=Qcold,phi=phi0,a=1.,alpha=0.,delay=0.)==min(S21_notch(f,fr,Ql=QLold,Qc=Qcold,phi=phi0,a=1.,alpha=0.,delay=0.)))[0,0]]
yhw=(1+S21min)/2

fend=f[np.argwhere(np.abs(np.abs(S21_notch_lhalf(f,a=1,f0=fr,Qe=2*Qcold,Qi=Qiold,phi=phi0,delay=0.))-yhw)==min(np.abs(np.abs(S21_notch_lhalf(f,a=1,f0=fr,Qe=2*Qcold,Qi=Qiold,phi=phi0,delay=0.))-yhw)))[0,0]]


plt.plot([fend-kappanew,fend],[yhw,yhw],"b-")
#plt.plot([fend-kappaold,fend],[yhw,yhw],"r--")
#plt.plot([fend-kappanew2,fend],[yhw,yhw],"k-")
