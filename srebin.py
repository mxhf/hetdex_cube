from scipy import interpolate
import numpy as np


def linlin(x, fluxv, start, step, stop=0., dowarn=True):
        """lin_rebin(x, flux, start, stop=0, step)
        Performs linear rebinning of arbitrarely sampled input spectra.
        x an array of input wavelengths in angstrom
        flux and array of fluxes in arbitrary untis
        start start wavelength in A
        stop (optional stop wavelength in A)
        step in km/s   
        function function for rebinning. Note for lin to log binning, "exp" has to be given, for 
                log to lin binning "log" has to be given!
        """
        
        warned = ~dowarn

        # make sure the array are sorted in ascending order     
        ii = np.argsort(x)
        x = x[ii]
        fluxv = fluxv[ii]
        
        if stop == 0. or stop > x[-1]:
                stop = x[-1]

        # prepare pix(wavelength) interpolation class   
        n1 = interpolate.interp1d( x, np.arange(x.shape[0]), 'linear', bounds_error=False)
        
        #output
        out = []
        
        # do the interpolation here, just one, this speeds the whole process up
        nmax = (stop - -start - step/2.)/step
        aalam  = np.arange(nmax+1)*step + start
        aalam1 = aalam + (-step/2.)
        aalam2 = aalam + (+step/2.)
        
        xx1  = n1(aalam1)
        xx2  = n1(aalam2)

        n = -1
        while True:
                # calculate current window size 
                n += 1
                #print "alam = %f + %f - %f " % (alam, alam1, alam2)
                #end while-loop if end of input range is reached
                #print(n)
                #print(aalam2)
                if aalam2[n] > stop:
                        break

                if aalam1[n] < x[0]:
                        out.append(0)
                        continue
                #tanslate wavelengths into pixel indices through linear interpolation
                x1  = xx1[n]
                x2  = xx2[n]
                mx1 = int( np.floor(x1) )
                mx2 = int( np.floor(x2) )
                if (x1 <= 0.) :
                        mx1 = 0
                if (x2 >= len(x)):
                        mx2 = len(x)-2
                # first integrate over all pixel which are fully contained in the
                #  current interval, integrate piecewise from the 
                #  center of on pixel to the next (therefore (pixel_i + pixel_(i+1))/2.)
                #
                flx=0.
                for ipix in range(int( mx1+1),int(mx2-1 + 1)):
                        flx=flx+( fluxv[ipix] + fluxv[ipix+1] )/2.
                
                # Now take care of the edges of the current interval using the
                # trapezium rule
                if (x1 < 0.):
                        x1 = 0.

                if (x2 > len(x)-1 ):
                         x2 = len(x)-1
                        
                if mx1 < mx2:
                        a1 = fluxv[mx1+1]-fluxv[mx1]
                        b1 = fluxv[mx1]
                        f_prime1 = a1*(x1-mx1) + b1
                        flx1 = (mx1+1-x1) * (f_prime1 + fluxv[mx1+1])/2.
                        
                        a2 = fluxv[mx2+1]-fluxv[mx2]
                        b2 = fluxv[mx2]
                        f_prime2 = a2*(x2 - mx2) + b2
                        flx2 = (x2-mx2) * ( f_prime2 + fluxv[mx2] )/2.
                elif mx1 == mx2:
                        if not warned:
                                print ("Warning: Such a small step width leads to sub pixel bins.")
                                warned = True                   
                        a1 = fluxv[mx1+1]-fluxv[mx1]
                        b1 = fluxv[mx1]
                        f_prime1 = a1*(x1-mx1) + b1
                        f_prime2 = a1*(x2-mx1) + b1
                        flx1 = (x2-x1) * (f_prime1+f_prime2)/2.
                        flx2 = 0.
                else:
                        print("Error mx1 > mx2, should never happen!")

                out.append(flx + flx1 + flx2)
        out_xs = np.arange(len(out))*step + (start)        

        return out_xs, np.array(out)

