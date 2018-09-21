library(deSolve)
library(FME)
library(manipulate)
library(readr)

data <- read.table("http://www.ndbc.noaa.gov/view_text_file.php?filename=ftpc1h2008.txt.gz&dir=data/historical/stdmet/")
Timeframe = data[60178:61967,]

#Cone parameters
a_sw = 0.86; #short-wave absorbivity 
a_lw = 0.90; #long-wave absorbivity
e_lws = 0.90; #long-wave emissivity of the shell
V_s = 0.7; #view factor
R = .0183; #Radius of shell
H = .0183; #Height of shell
K_r = 3.06; #Thermal conductivity of rock
p = 2601; #Density of rock
c_r = 789; #Specific heat of rock
k = 1.49*10^(-6); #Thermal diffusivity of rock
a = 1.66;  #
b = 0.389;
s = 5.67*10^(-8); #Stefan-Boltzmann constant

#Shortwave radiation  
T_a = Timeframe[,14] + 273.15#Temperature of the air(time dependent meteorlogical input) 
A_p =(1/4)*pi*R*sqrt(H^2+R^2);
x = seq(0, 8, length=length(T_a)) #Solar irradiance is estimated using triangle function (time dependent meteorological input)
t=0

trianglewave = 
  ifelse((x > 0 + t & x < (4/24) + t) | (x > 20/24 + t & x < 1 + t), 0,
         ifelse(x > (4/24) + t & x < (12/24) + t,  -((4/24 + t)-x),
                ifelse(x >= (12/24) + t & x < (20/24) + t, (20/24 + t)-x,0)))

plot(trianglewave)

t=1
while(t<9){
  trianglewave = trianglewave +
    ifelse((x > 0 + t & x < (4/24) + t) | (x > 20/24 + t & x < 1 + t), 0,
           ifelse(x > (4/24) + t & x < (12/24) + t,  -((4/24 + t)-x),
                  ifelse(x >= (12/24) + t & x < (20/24) + t, (20/24 + t)-x,0)))
  t=t+1 }

I_sw = trianglewave*(.9/trianglewave[112])
plot(I_sw)
q_1 =  A_p + a_sw + I_sw;

#Longwave radiation
A_l = pi*R*sqrt(H^2+R^2); 
e_lwa = (9.2*10^-7)*(T_a)^2 #long-wave emissivity of air for a clear sky
q_2 = V_s*A_l*e_lws*s*(T_a)^4*(e_lwa - 1); 
q_3 = 4*V_s*A_l*e_lws*s*(T_a)^3;

#Convective heat transfer  
A_cv = A_l;
K_a = 0.00501+(7.2*10^(-5))*T_a #conductivity of the air
v = (-1.25*10^(-5)) + (9.2*10^(-8))*T_a #kinematic viscosity of air
u = Timeframe[,7] #Wind speed(time dependent meteorological input)
u[125:126] = 1.5
h_c = a*((K_a*u/v)^b)*R^(b-1)#Convective heat transfer coefficient
q_4 = h_c*A_cv

#Conductive heat transfer
A_cd = pi*R^2
dx = .1;
q_5 = (A_cd*K_r)/(dx);
 
#T_2, the temperature a distance dx into the rock is estimated using a finite difference approach
L = 2
k = 1
dt = 8/length(T_a)
t = 8
T_o = Timeframe[,15] + 273.15 #Temperature of the ocean (time dependent meteoroloical input)
T_o[126] = 14.1 + 273.15

T_b = vector(length = t/dt)
T_2 = vector(length = t/dt)
T_b[1] = 288
T_2[1] = T_o[1]

m = vector(length = L/dx)
m[1] = T_b[1]
m[2:(L/dx)] = T_o[1]

Model = matrix(rep(0,len=length(m)*(t/dt)), nrow = L/dx)
Model[,1] = m
Model[1,] = m[1]
Model[(L/dx),] = m[(L/dx)]

for (j in 1:(length(T_b) - 1)){ 
  for (i in 2:(L/dx-1)){
    Model[i,j+1] = Model[i,j] + (k*dt/(dx)^2)*(Model[i+1,j]-2*Model[i,j]+Model[i-1,j])
  }
  T_2[j+1] = Model[2,j+1]
  T_b[j+1] = (q_1[j+1] + q_2[j+1] + (q_3[j+1]+q_4[j+1])*T_a[j+1] + q_5*T_2[j+1])/(q_3[j+1] + q_4[j+1] + q_5)
  Model[1,] = T_b[j+1]
  Model[(L/dx),] = T_o[j+1]
} 
plot(seq(0,8,length=length(T_b)),T_b)

plot(q_1)
plot(q_2)
plot(q_3)
plot(q_4)
plot(T_2)




