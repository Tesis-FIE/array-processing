#base simulation
import numpy as np
import matplotlib.pyplot as plt
rg = np.random.default_rng(1)
# Build a vector of 10000 normal deviates with variance 0.5^2 and mean 2
mu, sigma = 2, 0.5
noise = rg.normal(mu,sigma,10000)
plt.plot(noise)
plt.show()
N = 15
f = 2E9 #base frequency
c = 3E8 
wavelength = c/f  
n = 0:N-1 
n = n.' 
k0 = 2*pi / wavelength 
dtx = wavelength/2 

if (nargin == 1 )    # simple version
    startAngle = 0 
    endAngle = 60 
    N0 = 1   # Noise spectral density
else  
    N0 = 10^(-snr/10+log10(1)) 
end

 # # Defining angular space
dphi = linspace(0, 360, 10001) 
pas = zeros(1,length(dphi)) 

 # # Forming the continuous density
 # Broadside (90deg wrt array axis) is reference angle
 # Looking for the array index of a particular angle
p_i = find(abs(dphi -startAngle)<.018) 
pf = find(abs(dphi - endAngle)<.018) 
p_i = p_i(1) 
pf = pf(1) 
test = 1  
pas(p_i:pf) = test 
 

 # # Point signals
 #  test = find(abs(dphi - 80)<.1) 
 #  pas(test(1))= 1 
 #  second signal
 #  test = find(abs(dphi - 80)<.1) 
 #  pas(test(1))= 1 

 #Angular domain in radian
dphir = dphi.*(pi/180) 

 # Data acquisition
 # in this case, we must assume signals coming from the region described 
 # by nominal PAS

A = zeros(N, pf-p_i+1) 
A = exp(1j*k0*sin(dphir(p_i:pf))*dtx.*n)   # "continuous" paths or signals

 # We use only one snapshot thus x = A+Noise  t = 1 
 # The experiment is repeated Nr times
Nr = 50 
R = zeros(N) 
p = zeros(1,length(dphir)) 
pbb = zeros(1,length(dphir)) 

accP = zeros(1,length(dphir)) 
accPbb = zeros(1,length(dphir)) 

for i = 1:Nr
    Noise=randn(N,pf-p_i+1)*sqrt(N0) 
    x = A+Noise 
    R = x*x' 
    pNoise = sum(Noise.^2,2)./length(Noise) 
    pSignal = sum(abs(A).^2,2)./length(A) 
    snr = 10*log10(mean(pSignal)/mean(pNoise)) 
     # Then, an estimate of correlation matrix 
      #R = R./Nr 
     # Bartlett beamformer

    for i = 1:length(dphi)
        v= exp(1j*k0*sin(dphir(i))*dtx.*n) 
        pbb(i) = (v'*R*v)/(v'*v) 
    end
     # obtenido mediante SPAW
    for i = 1:length(dphi)
        v= exp(1j*k0*sin(dphir(i))*dtx.*n) 
        p(i) =k0*dtx/(2*pi)*abs(cos(dphir(i)))...
            * (v'*R*v)/(v'*v) 
    end
     #aplicado a mvdr
    
    accP = accP+p 
    accPbb = accPbb + pbb 
    
end
