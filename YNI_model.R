rm(list=ls())
library(FME)
library(deSolve)
library(gdata)
library(manipulate)


YNI = function(t,y,parms){
  with(as.list(c(y,parms)), {
    I_Na=c_Na*.5*m^3*h*(V-30)
    I_K=p*(.7*(exp(0.0277*(V+90))-1))/(exp(0.0277*(V+40)))
    I_l=0.8*(1-exp(-(V+60)/20))
    I_s=12.5*(0.95*d+0.05)*(0.95*f+0.05)*(exp((V-30)/15)-1)
    I_h=0.4*q*(V+25)
    alpha_m=(V+37)/(1-exp((V+37)/-10))
    beta_m=40*exp(-0.056*(V+62))
    alpha_h=1.209*10^-3*exp((V+20)/-6.534)
    beta_h=1/(1+exp((V+30)/-10))
    alpha_p=9*10^-3*1/(1+exp(-(V+3.8)/9.71))+6*10^-4
    beta_p=2.25*10^-4*(V+40)/(exp((V+40)/13.3)-1)
    alpha_q=3.4*10^-4*(V+100)/(exp((V+100)/4.4)-1)+4.95*10^-5
    beta_q=5*10^-4*(V+40)/(1-exp(-(V+40)/6))+8.45*10^-5
    alpha_d=1.045*10^-2*(V+35)/(1-exp(-(V+35)/2.5))+3.125*10^-2*V/(1-exp(-V/4.8))
    beta_d=4.21*10^-3*(V-5)/(exp((V-5)/2.5)-1)
    alpha_f=3.55*10^-4*(V+20)/(exp((V+20)/5.633)-1)
    beta_f=9.44*10^-4*(V+60)/(1+exp(-(V+29.5)/4.16))
    
    
    dV=(I_app-I_Na-I_K-I_l-I_s-I_h)*1/Cm
    dm=alpha_m*(1-m)-beta_m*m
    dh=alpha_h*(1-h)-beta_h*h
    dp=alpha_p*(1-p)-beta_p*p
    dd=alpha_d*(1-d)-beta_d*d
    df=alpha_f*(1-f)-beta_f*f
    dq=alpha_q*(1-q)-beta_q*q
    list(c(dV,dm,dh,dp,dd,df,dq));
  })
}
yini=c(V=-60,m=.10538, h=.82, p=0.301, d=2.25e-04, f=.795, q=.01022) # at V0=-40 no AP occurs
time=seq(0,1000,length=5000)
parms=c(I_app=.1, c_Na=1, Cm=1)
Soln=ode(y=yini,func=YNI,times=time, parms=parms)
#Soln
V=Soln[,2]
time=Soln[,1]
plot(time, V, type="l", main="Voltage vs. Time in YNI model", xlim=c(0,1000))

time=Soln[,1]
V=Soln[,2]
d=Soln[,6]
m=Soln[,3]
h=Soln[,4]
f=Soln[,7]
q=Soln[,8]
p=Soln[,5]


plot(time, V, type="l", main="Voltage vs. Time in YNI model", xlim=c(0,1000))

manipulate({
 parms=c(Cm=1, I_app=I, c_Na=1)
 Soln=ode(y=yini,func=YNI,times=time, parms=parms)
 plot(time, Soln[,2], type="l", main="YNI model Action Potentials", ylab="membrane potential (mV)",xlab="time (ms)", xlim=c(0,1000), ylim=c(-70,40))
},
I=slider(0,5,step=.05))

parms=c(Cm=1, I_app=0, c_Na=1)
Soln=ode(y=yini,func=YNI,times=time, parms=parms)
plot(time, Soln[,2], type="l", main="YNI model Action Potentials", ylab="membrane potential (mV)",xlab="time (ms)", xlim=c(0,1000), ylim=c(-70,40))
parms=c(Cm=1, I_app=1, c_Na=1)
Soln=ode(y=yini,func=YNI,times=time, parms=parms)
lines(time, Soln[,2],col=2)
#parms=c(Cm=1, I_app=1, c_Na=1)
#Soln=ode(y=yini,func=YNI,times=time, parms=parms)
#lines(time, Soln[,2],col=3)
parms=c(Cm=1, I_app=2, c_Na=1)
Soln=ode(y=yini,func=YNI,times=time, parms=parms)
lines(time, Soln[,2],col=3)
legend("topright", c("I=0","I=1","I=2"), col=c(1,2,3),lty=c(1,1,1))

#Soln[5000,]
c_Na=1

I_Na=c_Na*.5*m^3*h*(V-30)
I_K=p*(.7*(exp(0.0277*(V+90))-1))/(exp(0.0277*(V+40)))
I_l=0.8*(1-exp(-(V+60)/20))
I_s=12.5*(0.95*d+0.05)*(0.95*f+0.05)*(exp((V-30)/15)-1)
I_h=0.4*q*(V+25)
alpha_m=(V+37)/(1-exp((V+37)/-10))
beta_m=40*exp(-0.056*(V+62))
alpha_h=1.209*10^-3*exp((V+20)/-6.534)
beta_h=1/(1+exp((V+30)/-10))
alpha_p=9*10^-3*1/(1+exp(-(V+3.8)/9.71))+6*10^-4
beta_p=2.25*10^-4*(V+40)/(exp((V+40)/13.3)-1)
alpha_q=3.4*10^-4*(V+100)/(exp((V+100)/4.4)-1)+4.95*10^-5
beta_q=5*10^-4*(V+40)/(1-exp(-(V+40)/6))+8.45*10^-5
alpha_d=1.045*10^-2*(V+35)/(1-exp(-(V+35)/2.5))+3.125*10^-2*V/(1-exp(-V/4.8))
beta_d=4.21*10^-3*(V-5)/(exp((V-5)/2.5)-1)
alpha_f=3.55*10^-4*(V+20)/(exp((V+20)/5.633)-1)
beta_f=9.44*10^-4*(V+60)/(1+exp(-(V+29.5)/4.16))

plot(time, V, type="l", main="Voltage vs. Time in YNI model", xlim=c(0,450), xlab="time (ms)")
lines(time, I_K, col=2)
lines(time, I_l, col=3)
lines(time, I_s, col=4)
lines(time, I_h, col=5)


plot(time, I_h, ylab="I_ion (microAmps/cm^2)", xlab="time (ms)", type="l", ylim=c(-6,2), xlim=c(20,600),main="Ion Currents versus time")
lines(time, I_K, col=2)
lines(time, I_l, col=3)
lines(time,I_s, col=4)
lines(time, I_Na, col=6)
legend("bottomright",c("I_h","I_K","I_l","I_s","I_Na"),col=c(1,2,3,4,6),lty=c(1,1,1,1,1))

tau_h=1/(alpha_h+beta_h)
tau_m=1/(alpha_m+beta_m)
tau_p=1/(alpha_p+beta_p)
tau_q=1/(alpha_q+beta_q)
tau_d=1/(alpha_d+beta_d)
tau_f=1/(alpha_f+beta_f)
plot(time, tau_f, ylab="time constants tau",col=4,main="gate time constants versus time",type="l", ylim=c(-150,400), xlim=c(0,320))
lines(time, tau_m, col=2)
lines(time, tau_p, col=6)
lines(time, tau_q, col=5)
lines(time, tau_d, col=1)
lines(time, tau_h, col=3)
legend("topleft",c("f","m","p","q","d","h"),col=c(4,2,6,5,1,3),lty=c(1,1,1,1,1,1))

plot(time, d, main="Ion gates versus time", ylab="gates", type="l", xlim=c(0,350), xlab="time (ms)")
lines(time, m, col=2)
lines(time, h, col=3)
lines(time, f, col=4)
lines(time, q, col=5)
lines(time, p, col=6)
legend("topleft", c("d","m","h","f","q","p"), col=c(1, 2, 3, 4, 5, 6), lty=c(1,1,1,1,1,1))

plot(V, d, main="Ion gates versus membrane potential",type="l", ylab="gating variables", xlab="membrane potential (mV)")
lines(V, m, col=2)
lines(V, h, col=3)
lines(V, f, col=4)
lines(V, q, col=5)
lines(V, p, col=6)
legend("topright", c("d","m","h","f","q","p"), col=c(1, 2, 3, 4, 5, 6), lty=c(1,1,1,1,1,1))

plot(V, tau_f, col=4,xlab="membrane potential (mV)", ylab="tau", main="time constants versus membrane potential", type="l", ylim=c(0,350))
lines(V, tau_m, col=2)
lines(V, tau_p, col=6)
lines(V, tau_q, col=5)
lines(V, tau_d, col=1)
lines(V, tau_h, col=3)
legend("topleft",c("f","m","p","q","d","h"),col=c(4,2,6,5,1,3),lty=c(1,1,1,1,1,1))

# Steady state functions

d_inf=alpha_d/(alpha_d+beta_d)
m_inf=alpha_m/(alpha_m+beta_m)
h_inf=alpha_h/(alpha_h+beta_h)
f_inf=alpha_f/(alpha_f+beta_f)
q_inf=alpha_q/(alpha_q+beta_q)
p_inf=alpha_p/(alpha_p+beta_p)
plot(V, d_inf, ylab="equilibrium gate values", main="Steady State gate values",type="l",xlab="membrane potential")
lines(V, m_inf, col=2)
lines(V, h_inf, col=3)
lines(V, f_inf, col=4)
lines(V, q_inf, col=5)
lines(V, p_inf, col=6)
legend("right", c("d_inf","m_inf","h_inf","f_inf","q_inf","p_inf"), col=c(1, 2, 3, 4, 5, 6), lty=c(1,1,1,1,1,1))

#Bifurcation of I
manipulate({
  c_Na=1
  Cm=1
  I_Na=c_Na*.5*m^3*h*(V-30)
  I_K=p*(.7*(exp(0.0277*(V+90))-1))/(exp(0.0277*(V+40)))
  I_l=0.8*(1-exp(-(V+60)/20))
  I_s=12.5*(0.95*d+0.05)*(0.95*f+0.05)*(exp((V-30)/15)-1)
  I_h=0.4*q*(V+25)
  alpha_m=(V+37)/(1-exp((V+37)/-10))
  beta_m=40*exp(-0.056*(V+62))
  alpha_h=1.209*10^-3*exp((V+20)/-6.534)
  beta_h=1/(1+exp((V+30)/-10))
  alpha_p=9*10^-3*1/(1+exp(-(V+3.8)/9.71))+6*10^-4
  beta_p=2.25*10^-4*(V+40)/(exp((V+40)/13.3)-1)
  alpha_q=3.4*10^-4*(V+100)/(exp((V+100)/4.4)-1)+4.95*10^-5
  beta_q=5*10^-4*(V+40)/(1-exp(-(V+40)/6))+8.45*10^-5
  alpha_d=1.045*10^-2*(V+35)/(1-exp(-(V+35)/2.5))+3.125*10^-2*V/(1-exp(-V/4.8))
  beta_d=4.21*10^-3*(V-5)/(exp((V-5)/2.5)-1)
  alpha_f=3.55*10^-4*(V+20)/(exp((V+20)/5.633)-1)
  beta_f=9.44*10^-4*(V+60)/(1+exp(-(V+29.5)/4.16))
  V=seq(-90,40,length=5000);
  dV=(I-I_Na-I_K-I_l-I_s-I_h)*1/Cm
  plot(V, dV, type="l",ylab="dV/dt")
  abline(0,0)
},
I=slider(0,5,step=.05))

Vvals=seq(-120,50,length=5000)
I_ap=(I_Na-I_K-I_l-I_s-I_h)
plot(I_ap, Vvals, type="l", xlim=c(-1,1), ylim=c(-30,10),main="Membrane potential equilibrium points (YNI model)", ylab="membrane potential (mV)", xlab="Applied current")

