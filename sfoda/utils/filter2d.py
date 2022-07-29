import numpy as np
from scipy.fftpack import fft2, ifft2, fftfreq, fftshift, ifftshift

###
# Filter helper functions
def butterworth(f, f0, n):
    omega = f/f0
    return np.power(np.sqrt(1+omega**(2*n)),-1)

def butterworth_highpass(f, f0, n):
    
    return 1 - butterworth(f, f0, n)

def butterworth_bandpass(f, f1, f2, n):
    
    H1 = butterworth(f, f1, n)
    H2 = butterworth_highpass(f, f2, n)
    return H1*H2

###
# Main function
def filter2d(Z3, dx, l_low, l_high, thetalow, thetahigh, n=8):
    """
    Two-dimensional wavenumber-directional Fourier transform filter.
    Uses a Butterworth filter.
    
    Inputs:
    ------
        z: 2D complex array
        dx: grid spacing
        l_low: shortest wavelength
        l_high: longest wavelength
        thetalow: low angle for filter (degrees CCW from E)
        thetahigh: high angle for filter (degrees CCW from E)
        
        n: (optional, default=8) order of the Butterworth filter
        
    Outputs:
    -----
        zf: 2D filtered complex array
        H: the filter function
        kx, ky: horizontal wavenumbers
    """        
    Z = fft2(Z3)

    ky = fftfreq(Z3.shape[0], dx/(2*np.pi))
    kx = fftfreq(Z3.shape[1], dx/(2*np.pi))
    kx = fftshift(kx)
    ky = fftshift(ky)
    
    Lx,Ly = np.meshgrid(kx,ky)
    K = np.abs(Lx + 1j*Ly)

    # Wavenumber filter
    K_low = 2*np.pi/l_low
    K_high = 2*np.pi/l_high

    Hk = butterworth_bandpass(K, K_low, K_high, n)
     
    degrad = np.pi/180
    theta = np.angle(Lx + 1j*Ly)
    thetadeg = np.mod(theta*180/np.pi,360)
    theta = np.mod(theta, 360*degrad)

    ###
    # Create the filter matrix
    
    # Rotate everything so thetamid is 180
    thetamid = 0.5*(thetalow+thetahigh)
    thetarotate = 180-thetamid

    Hf1 = butterworth_highpass(theta+thetarotate*degrad, thetalow*degrad+thetarotate*degrad, n)
    Hf2 = butterworth(theta+thetarotate*degrad, thetahigh*degrad+thetarotate*degrad,n)
    Hf = Hf1*Hf2
    
    # Combine the filters
    H =  Hk*Hf
    
    # Now reorder H into the original FFT ordering
    H_r = ifftshift(H,axes=1)
    H_r = ifftshift(H_r,axes=0)

    # Finally, filter
    Zf = ifft2(Z*H_r)
    
    return Zf, H, kx, ky