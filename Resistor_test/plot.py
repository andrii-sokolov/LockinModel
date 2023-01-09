import matplotlib.pyplot as plt
import numpy as np

[f1,X1,Y1,R1,Theta1] = np.genfromtxt("test_1.csv",skip_header=1,delimiter=",").T
GTiA1 = 1.0e-9 #[A/V]
Vin_rms1 = 1e-3*np.sqrt(2) #V
[f2,X2,Y2,R2,Theta2] = np.genfromtxt("test_2.csv",skip_header=1,delimiter=",").T
GTiA2 = 1.0e-9 #[A/V]
Vin_rms2 = 1e-3*np.sqrt(2)
[f3,X3,Y3,R3,Theta3] = np.genfromtxt("test_3.csv",skip_header=1,delimiter=",").T
GTiA3 = 1.0e-9 #[A/V]
Vin_rms3 = 1e-3*np.sqrt(2)
[f4,X4,Y4,R4,Theta4] = np.genfromtxt("test_4.csv",skip_header=1,delimiter=",").T
GTiA4 = 1.0e-9 #[A/V]
Vin_rms4 = 1e-3*np.sqrt(2)

fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f1,R1*GTiA1/Vin_rms1*1e12,".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f1,Theta1,"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.plot([f1[0],f1[-1]],[1,1],"--",linewidth=0.5)
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$\\Theta$, $^0$")
axs1.set_ylabel("$G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("test1.pdf")
plt.show()

fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f2,R2*GTiA2/Vin_rms2*1e12,".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f2,Theta2,"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.plot([f1[0],f1[-1]],[1,1],"--",linewidth=0.5)
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$\\Theta$, $^0$")
axs1.set_ylabel("$G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("test2.pdf")
plt.show()

fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f3,R3*GTiA3/Vin_rms3*1e12,".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f3,Theta3,"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.plot([f3[0],f3[-1]],[1,1],"--",linewidth=0.5)
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$\\Theta$, $^0$")
axs1.set_ylabel("$G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("test3.pdf")
plt.show()

fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f4,R4*GTiA4/Vin_rms4*1e12,".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f4,Theta4,"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.plot([f4[0],f4[-1]],[1,1],"--",linewidth=0.5)
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$\\Theta$, $^0$")
axs1.set_ylabel("$G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("test4.pdf")
plt.show()


fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f1,np.abs(R1*GTiA1/Vin_rms1*1e12-1),".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f1,np.abs(Theta1),"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$|\\Theta|$, $^0$")
axs1.set_ylabel("$\\Delta G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("err1.pdf")
plt.show()

fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f2,np.abs(R2*GTiA2/Vin_rms2*1e12-1),".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f2,np.abs(Theta2),"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$|\\Theta|$, $^0$")
axs1.set_ylabel("$\\Delta G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("err2.pdf")
plt.show()

fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f3,np.abs(R3*GTiA3/Vin_rms3*1e12-1),".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f3,np.abs(Theta3),"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$|\\Theta|$, $^0$")
axs1.set_ylabel("$\\Delta G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("err3.pdf")
plt.show()

fig1, axs1 = plt.subplots()
x2 = axs1.twinx()
x2.set_xscale("log")
p1 = axs1.plot(f4,np.abs(R4*GTiA4/Vin_rms4*1e12-1),".-",linewidth=0.5, label="$\\leftarrow$ $G$")
p2 = x2.plot(f4,np.abs(Theta4),"r--",linewidth=0.5, label="$\\rightarrow$ $\\Theta$")
axs1.set_xlabel("$f$, Hz")
x2.set_ylabel("$|\\Theta|$, $^0$")
axs1.set_ylabel("$\\Delta G$, pSi")
lns = p1+p2
labs = [l.get_label() for l in lns]
axs1.legend(lns, labs, loc=0)
axs1.grid()
plt.savefig("err4.pdf")
plt.show()

from scipy import signal

err1 = np.abs(R1*GTiA1/Vin_rms1*1e12-1)
err2 = np.abs(R2*GTiA2/Vin_rms2*1e12-1)
err3 = np.abs(R3*GTiA3/Vin_rms3*1e12-1)
err4 = np.abs(R4*GTiA4/Vin_rms4*1e12-1)

corr1 = signal.correlate((err1 - np.mean(err1)) / (np.std(err1)* len(err1)), (np.abs(Theta1) - np.mean(np.abs(Theta1))) / (np.std(np.abs(Theta1))), mode='same')
corr2 = signal.correlate((err2 - np.mean(err2)) / (np.std(err2)* len(err2)), (np.abs(Theta2) - np.mean(np.abs(Theta2))) / (np.std(np.abs(Theta2))), mode='same')
corr3 = signal.correlate((err3 - np.mean(err3)) / (np.std(err3)* len(err3)), (np.abs(Theta3) - np.mean(np.abs(Theta3))) / (np.std(np.abs(Theta3))), mode='same')
corr4 = signal.correlate((err4 - np.mean(err4)) / (np.std(err4)* len(err4)), (np.abs(Theta4) - np.mean(np.abs(Theta4))) / (np.std(np.abs(Theta4))), mode='same')

plt.plot(f1,corr1,label="Experiment 1")
plt.plot(f2,corr2,label="Experiment 2")
plt.plot(f3,corr3,label="Experiment 3")
plt.plot(f4,corr4,label="Experiment 4")
plt.xscale("log")
plt.xlabel("$f$, Hz")
plt.ylabel("Correlation function")
plt.legend()
plt.grid()
plt.savefig("correlation.pdf")
