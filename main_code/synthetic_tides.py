# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 16:13:30 2015

@author: Basile HECTOR

Module for tides calculations functions (solid-earth and ocean loading).
Adapted from Matlab(TM) codes by R. Cattin:

Cattin, R., Mazzotti, S., and Baratin, L.-M., 2015, 
GravProcess: An easy-to-use MATLAB software to process campaign gravity data 
and evaluate the associated uncertainties: Computers & Geosciences,
 v. 81, p. 20–27, doi: 10.1016/j.cageo.2015.04.005.


to be modified or checked:
- global variables or class attribute?
- altitude = 0: remove the hard-wire
- beware, in the Matlab code, gl is never assigned which makes the code to fail
  in the detail, we find somewhere g1=g; which should (certainly) be gl 
"""
import numpy as np
from datetime import *
from matplotlib.dates import date2num,num2date
import math

def earth_tide(theta,lamda,gtime):
    """
    [TIDE] = earth_tide(LAT,LON,gtime)
    input
    LAT,LON - North Latitude and East Longitude, degrees
    gtime - a single time in the datetime python format
    output:
    computes the tidal correction in microgals.
    it should be added to your measurements. 
    
    
    GDC returns the constant tide (Honkasalo term) as well.  The
    constant tide is included in TIDE. 
    In order to correct surface gravity records for tides, ADD this 
    tide output to gravity records
    
    BOMM subroutine written by Jon Berger  November 1969
    astronomy revised by Judah Levine (after J C Harrison) Sept 1973
    updated by Karen Young  March 1977, for PDP 11
    amended by Duncan Agnew June 1978 to speed execution
    solar astronomy redone and rigid-earth potential added by Duncan
    Agnew June 1979
    tide generating part put in separate subroutine and computation
    of Munk-Cartwright coefficients added by Duncan Agnew Feb 1981,
    uly 1982
    This version rewritten for export, using F77 calls for I/O,
    by Duncan Agnew Mar 1987
    
    This version stripped down to theoretical gravity tide and
    called as a subroutine.  Output is in a passed vector.
    By Glenn Sasagawa, June 1988 
     
    Code modified to Matlab(tm) by R. Cattin in 2015

    The present version has been translated from Matlab to Python 
    by B. Hector in 2015.    
    
    tides are calculated from harmonics 2 through 3 for
    the lunar terms, 2 for the solar terms.
    love numbers h(n), k(n), and l(n) are set by data statements.
     
    gravity tide is in microgals, plus for up acceleration
     
    Arguments
    yr1,day1,zhr	Start year, day and hour
    yr2,day2,yhr	Start year, day and hour
    NOTE THAT ZHR AND YHR ARE DOUBLE PRECISION
    d			Time interval in hours, DOUBLE PRECISION
    theta		North latitude
    lamda		East longitude
    iterms		Number of output tide terms
    gravtide		Output gravity tide    
    
    """
        
        
    global dsz, dcz, dsl, dcl, ssz, scz, ssl, scl, dpar, sdist  # bpos common block
    global h, k, l                	# love common block
    h =[0.6114,0.2891,0.175]
    k =[0.304,0.09421,0.043]
    l =[0.0832,0.0145,0.0103]
    
    global azt, azs        # azimut common block
    global etmut                  # tdiff common block
    global moon                  # sunny common block
    moon=0       
    #hardwire these - you can only send it ONE droptime
    deltat = 1
    NPT = 1
    YY = gtime.year
    MO = gtime.month
    DD = gtime.day
    HH = gtime.hour
    MM = gtime.minute
    SS = gtime.second
    # Initialize variables          
    irl=1
    iflag=0
    ntotl=1
    iget=[0, 0, 0, 0, 0, 0, 0] # ' !!!
    ispc=[0, 0, 0, 0]# ' !!!
    ntw=[1, 0, 0]# ' !!!
    ioptn ='t'
    ielement = 0 
    #	data statements for input and output unit numbers (on terminal I/O)
    inun=5
    ioun=6
    nptpb=6

#       
    yr1=YY-1900
    day1 = date2num(datetime(YY,MO,DD))
    #	find times in hours from 0 hr, 1 jan 1900
    #matlab:
    ts=SS/3600+MM/60 + HH+24*(day1-1)+8760*yr1+24*np.fix((yr1-1)/4)
    # python:
    dj=date_to_julian_day(datetime(YY,MO,DD))
    djref=date_to_julian_day(datetime(1899,12,31,0,0,0))
    delta_dj=dj-djref #difference in days from current date (0hr) to 0hr, 1 jan 1900
    delta_djhr=float(delta_dj)*24.0+HH-12.0+MM/60.0+SS/3600.0
    te = ts + (NPT-1)*deltat/3600
    d = deltat/3600
    #terms=(te-ts)/d + 1 
    terms=NPT
     
    #done asking questions - begin execution
    i = 1
    tt = ts
    sph(theta,lamda,0)
    etmut = 41.184 + yr1 - 70
    #matlab:
    #t = (tt+12 + (etmut/3600))/876600
    t = (delta_djhr+ etmut/3600)/876600
    #  t is ephemeris time in julian centuries from 12 hr 0 jan 1900
    ephem(t)
    
    # calculate normalized gravity tides
    [grav, tilt, strain, gdc] = elastd(ntw)

    gravtide=1.e5*grav
    #convert m/s² to mgal: 1m/s² = 100 gal = 100 000 mgal
    
    iflag = 1
    
    iterms = np.fix(terms)
    i = 1
    return gravtide

def sph(grlat,elong,ht):
    """    
    #   for a point at geographical north latitude grlat, east longitude elong      
    #   (in degrees), and height ht (in meters), finds geocentric position   
    #   and local g using formulae for a spheroid
    """
    
    # Initialize Variables 
    global cth, sth, clg, slg, dif, radn, gl 		#common/obs/
    gn=9.798277692
    ae=6378140.
    f=0.00335281
    rm=0.00344978
    dr=0.01745329252    
    
    clong = np.cos(elong*dr)    
    slong = np.sin(elong*dr)    
    # latitude difference 
    dvert = f*(1.+.5*f)*np.sin(2.*grlat*dr) - .5*f*f*np.sin(4.*grlat*dr)     
    gcclat = (3.1415926535898/2.) - (grlat*dr - dvert)   
    cthet = np.cos(gcclat)      
    sthet = np.sin(gcclat)      
    # geocentric radius   
    radn = 1 - f*(cthet**2)*(1 + 1.5*f*(sthet**2))     
    # formulae for g are from jeffreys, 4.022 and 4.023      
    g = gn*(1 + f - 1.5*rm + f*(f-(27/14)*rm) + (2.5*rm - f -f*(f-\
          (39/14)*rm))*(cthet**2) - (f/2)*(7*f-15.*rm)*((cthet*sthet)**2))      
    # free air correction 
    g = g - g*(2.*ht*(1.+f+rm - 2.*f*(cthet**2))/ae)     
    
    # Conversion Here for Globals
    cth=cthet
    sth=sthet
    clg=clong
    slg=slong
    dif=dvert
    gl=g
    
    
    
def ephem(t):
    """
    t is ephemeris time in julian centuries from 12 hr 0 jan 1900
    (in ordinary civil reckoning this is greenwich noon on 31 december
    1899). if the difference between ephemeris and unversal time is
    not put in (see common tdiff below), t should be in universal time
    which is (nearly) the time ordinarily kept by clocks.
    computes positions of the sun and moon at time t, returning results
    in common block bpos. the solar ephemeris uses the mean sun.
    Derived from J. Levine's revision (after J. C. Harrison)
    of an earthtide program by J. Berger and W. E. farrell, with small
    alterations by D. C. Agnew, partly after M. Wimbush. present
    subroutine version by d. c. agnew.
     
    common block bpos contains, in order:
    sine and cosine of colatitude of sublunar point
    sine and cosine of east longitude of sublunar point
    sine and cosine of colatitude of subsolar point
    sine and cosine of east longitude of subsolar point
    the lunar sine parallax in arc seconds
    the solar distance in astronomical units
    """
    global dsz, dcz, dsl, dcl, ssz, scz, ssl, scl, dpar, sdist # bpos common block
    #
    #  common block containing the difference ephemeris minus
    #  universal time, in seconds. if this is not known it should
    #  be set to zero, and the argument to the program should be
    #  universal rather than ephemeris time.
    #
    global etmut #tdiff common block
    #  common block containing the instruction on which ephemeris to compute
    #  moon =   0  - both sun and moon
    #      1  - moon only
    #      2  - sun only
    global moon

    pi20=62.8318530717958
    #   compute universal time in hours
    ts = 876600*t -12 - (etmut/3600)
    hr = np.mod(ts,24)
    #   compute obliquity of the ecliptic
    w = .409319747 - .0002271107*t
    cosw = np.cos(w)
    sinw = np.sin(w)
    t2 = t*t
    if (moon!=1):
        # compute solar constants for given t
        hs = 4.881627982482 + 628.3319508731*t + 0.523598775578*10**(-5)*t2
        hs = np.mod(np.mod(hs,pi20)+pi20,pi20)
        ps = 4.908229466993+ 0.03000526416690*t + 0.790246300201*10**(-5)*t2
        es = 0.01675104 - 0.00004180*t - 0.000000126*t2
        psig = 0.2617993877971*(hr-12.) + hs
        chmp = np.cos(hs-ps)
        shmp = np.sin(hs-ps)
        ls = hs + shmp*es*(2.+2.5*es*chmp)
        sls = np.sin(ls)
        cz = sinw*sls
        sz = np.sqrt(1.-cz**2)
        psis = math.atan2(cosw*sls,np.cos(ls))
        rbarr = 1. + es*(chmp + es*(chmp-shmp)*(chmp+shmp))
        ll = psis - psig
        scz = cz
        ssz = sz
        ssl = np.sin(ll)
        scl = np.cos(ll)
        sdist = 1/rbarr
    
    # compute lunar constants for given t
    
    if (moon==2): 
       return 
    hm=4.7199666+8399.7091449*t-.0000198*t2
    pm=5.83515154+71.01804120839*t-.180205*10**(-3)*t2
    nm=4.523601515-33.75714624*t+.3626406335*10**(-4)*t2
    #   bl bls bf bd are the fundamental arguments of browns theory
    bl=hm-pm
    bls=hs-ps
    bf=hm-nm
    bd=hm-hs
    #   lunar lat long and parallax from brown.  latter two from
    #   improved lunar ephemeris, latitude from ras paper of 1908...
    tlongm=hm+.10976*np.sin(bl)-.02224*np.sin(bl-2.*bd)+0.01149*np.sin(2.*bd)+\
           0.00373*np.sin(2.*bl)-.00324*np.sin(bls)-.00200*np.sin(2.*bf)-0.00103*\
           np.sin(2.*bl-2.*bd)-.00100*np.sin(bl+bls-2.*bd)+0.00093*np.sin(bl+2.*bd)\
           -.00080*np.sin(bls-2.*bd)+.00072*np.sin(bl-bls)-.00061*np.sin(bd)-\
           .00053*np.sin(bl+bls)
    tlatm=.08950*np.sin(bf)+.00490*np.sin(bl+bf)-.00485*np.sin(bf-bl)-.00303*\
          np.sin(bf-2.*bd)+.00097*np.sin(2.*bd+bf-bl)-.00081*np.sin(bl+bf-2.*bd)+\
          .00057*np.sin(bf+2.*bd)
    plx=(3422.45+186.54*np.cos(bl)+34.31*np.cos(bl-2.*bd)+28.23*np.cos(2.*bd\
         )+10.17*np.cos(2.*bl)+3.09*np.cos(bl+2.*bd)+1.92*np.cos(bls-2.*bd)+1.44*\
         np.cos(bl+bls-2.*bd)+1.15*np.cos(bl-bls)-0.98*np.cos(bd)-0.95*np.cos(bl+\
         bls)-0.71*np.cos(bl-2.*bf)+0.62*np.cos(3.*bl)+0.60*np.cos(bl-4.*bd))
    sinmla=np.sin(tlatm)
    cosmla=np.cos(tlatm)
    sinmln=np.sin(tlongm)
    cosmln=np.cos(tlongm)
    #...convert from celestial lat and long according to explan suppl of
    #......na and le page 26
    cz=cosmla*sinmln*sinw+sinmla*cosw
    sz=np.sqrt(1.-cz**2)
    at1=cosmla*sinmln*cosw-sinmla*sinw
    at2=cosmla*cosmln
    ram=math.atan2(at1,at2)
    ll=ram-psig
    dcz = cz
    dsz = sz
    dsl = np.sin(ll)
    dcl = np.cos(ll)
    dpar = plx
          
    #------------------------------------------------------------------      
def elastd(ntw):
    """
    #   computes the earth tides on an elastic earth, given the solar and lunar
    #   positions and the place of observation.  degree 2 is used for the solar
    #   tides, 2 through 4 for the lunar tides. the results are returned in
    #   grav, tilt, and strain. ntw(1) gives the number of gravity tides
    #   ntw(2) the number of tilt tides, ntw(3) the number of strain tides.
    #   as the dimensioning shows, at most one gravity tide, two tilts, and
    #   three strains may be computed.  if ntw(1) = -1, the program will
    #   put the equilibrium potential height for a rigid spherical earth in
    #   grav.  the units are m/s**2, radians, extension, and meters.
    #   the sign convention is positive for upward potential
    #   height, a decrease in g, a tilt in the azimuth given, and
    #   extensional strain.
    #   that part of the tidal signal which comes from the permanent
    #   deformation is subtracted out using the coefficient of .31455 m
    #   found by Cartwright and Edden for the dc tide.
    #
    #   based (very closely) on J. Berger earth tide program.
    #   converted to subroutine by D. Agnew Nov 1979.
    #   computation of the potential and subtraction of dc tides added
    #   by D. Agnew Jan 1980.
    #
    """
    # Simulated common block obs with observatory information
    global cth, sth, clg, slg, dif, radn, gl 
    # Simulated common block love with love numbers 
    global  h, k, l      
    # Simulated common block azimut with strainmeter and tiltmeter azimuths
    global azt, azs 
    # Simulated common block bpos with lunar and solar colat and long, lunar sine parallax,
    # and solar distance
    global dsz, dcz, dsl, dcl, ssz, scz, ssl, scl, dpar, sdist  
    
    coor =[dsz, dcz, dsl, dcl, ssz, scz, ssl, scl]
    par = dpar 	
    # Data for mean parallaxes, a times m over m(earth), equatorial radius.
    rbor=[1.6592496e-2, 4.2635233e-5]
    amrat=[78451.25, 2.1235762e12]
    a=6.37814e6
    g1=9.79828
    g2=9.82022 
    ppp=[3,0,0]
    iflag=0 
    strain = [0]
    dele=[0,0,0]
    dim=[0,0,0]
    # On first call compute factors for gravity and tilt, and dc tides
    # at the given latitude.
    if (iflag != 1):
       iflag=1
       for i in range(3):
           #browse 0,1,2
           dele[i] = 1. + (2./(i+2.))*h[i] - ((i+3.)/(i+2.))*k[i]
           dim[i] = 1. + k[i] - h[i]
        #  dc gravity tide is also known as the Honkasalo correction
        #  **note that the love numbers for an elastic earth are used
        #  in computing the dc tide as well.eq
       gdc = -3.0481e-7*(3*cth**2 - 1.)*dele[0]*radn
       tnsdc = -9.1445e-7*cth*sth*dim[0]*radn/gl
       etdc = -1.555e-8*(h[0]*(3.*cth**2-1.) - 6.*l[0]*(2.*cth**2 -1.))
       eldc = -1.555e-8*(h[0]*(3.*cth**2-1.) - 6.*l[0]*cth**2)
       potdc = .0992064*(1.-3*cth**2)
       re = 1./(radn*a)
    
    # zero out arrays
    tilt = [0, 0]
    e = [0, 0, 0]
    tltcor = [0, 0]
    grav = 0
    gnth = 0
    # compute normalized parallax
    pa=[par/3422.45,1/sdist]

    # in outer loop, ii = 0 for moon, 1 for sun (basile)
    for ii in [0,1]:
        id = 3
        if(ii==1):
            id = 1 
        ir = 4*(ii)
        # find cosine of zenith angle, potential constants, legendre polynomials
        # and their derivatives, and derivatives of the cosine of the zenith angle.
        cll = clg*coor[ir+3] + slg*coor[ir+2]
        sll = slg*coor[ir+3] - clg*coor[ir+2]
        cz = coor[ir+1]
        sz = coor[ir]
        cu = cth*cz + sth*sz*cll
        xi = rbor[ii]*pa[ii]*radn
        cc = amrat[ii]*rbor[ii]*pa[ii]
        rkr=[cc*xi*xi,0,0]
        rkr[1] = rkr[0]*xi
        rkr[2] = rkr[1]*xi        

        p=[0.5*(3*cu*cu - 1.),0,0]
        pp=[3*cu,0,0]
        if(ii!= 1):
            p[1] = .5*cu*(5.*cu*cu - 3.)
            p[2] = .25*(7.*cu*p[1] - 3.*p[0])
            pp[1] = 1.5*(5.*cu*cu - 1.)
            pp[2] = .25*(7.*p[1] + cu*pp[1]) - 3.*pp[0]
            ppp[1] = 15.*cu
            ppp[2] = 7.5*(7.*cu*cu - 1.)        
        

        cut = -sth*cz + cth*sz*cll
        cutt = -cu
        cul = -sth*sz*sll
        cull = -sth*sz*cll
        cutl = -cth*sz*sll
        #for j = 1:id:
        for j in range(id):
            if(ntw[0]==1):
                grav = grav + dele[j]*(j+2)*rkr[j]*p[j]*g1*re
        gnth = gnth - dim[0]*rkr[0]*pp[0]*g1*cut*re
    # ellipticity corrections, convert strains to strainmeter
    if (ntw[0]==1):
        grav = grav + gnth*dif - gdc
    
    return [grav,tilt,strain,gdc]     
    
def date_to_julian_day(my_date):
    """
    Returns the Julian day number of a date.
    http://code-highlights.blogspot.fr/2013/01/julian-date-in-python.html
    number of days since November 24, 4714 BC
    """
    a = (14 - my_date.month)//12
    y = my_date.year + 4800 - a
    m = my_date.month + 12*a - 3
    return my_date.day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045

def oceanLoading(time,amp,phases,lon):
    """
    time is in hours from Jan 1, 1900, stationlongitude is longitude
    
        purpose:  To assess the ocean loading/gravity effect due to the
                  separate tidal primary components derived from the Agnew
                  program ( LOADF ).
    produce the ocean loading correction in µgal (?)
    
    Arguments:
    - time                 current time (datetime object)
    - amplitudes & phases: lists with amplitudes and phases. List order is:
        M2,S2,K1,O1,N2,P1,K2,Q1,Mf,Mm,Ssa      
        (amp,phases)
    - lon:     site longitude                          
       
    From R. Cattin, translated from MATLAB(TM) to Python by B. Hector
    """
    #ctime should be current time (in hours from Jan 1, 1900)
    djref=date_to_julian_day(datetime(1900,1,1,0,0,0))
    dj=date_to_julian_day(datetime(time.year,time.month,time.day))
    ctime=float(dj-djref)*24.0+time.hour+time.minute/60.0+time.second/3600.0

    dtor=np.pi/180
    TL=64.184       # time shift ET-UT = 41.184+yr-1970 seconds
    ets=(ctime/24+0.5+TL/86400)/36525-1;
    #Calculate the mean longitudes of the sun (h0),
    #moon (s0) and lunar perigee (p0).
    #And local siderial time (tl).
    #All quantities in radians
    p0 = lunarperigee(ets)
    h0 = sunlongitude(ets)
    s0 = moonlongitude(ets)
    AM = 18.69735+2400.05130*(ets-TL/3.15576e9)
    P = 0.00029*np.cos((1934*ets+145)*dtor)
    gast = np.pi+(ctime+AM+P)*np.pi/12 + lon*dtor
    return gravityeffect(h0,s0,p0,gast,amp,phases)
    
def lunarperigee(time):
    """
    LONGITUDE OF LUNAR PERIGEE IN RADIANS    
    TIME:  (JED-2451545.0)/36525             EPOCH = J2000.0
    ets is the Julian centuries since Jan 1.5 2000.
    t1 is the Julian centuries since Jan 0.5 1900.
    """
    dtor=np.pi/180
    t1=1+time
    t2=t1*t1
    t3=t2*t1
    perigee= 334.329653*dtor + 4069.0340329575*dtor*t1 - 0.010325*dtor*t2 - 1.2e-5*dtor*t3
    return perigee
    
def sunlongitude(time):
    """    
    LONGITUDE OF THE SUN IN RADIANS ( NUTATION IS TAKEN INTO ACCOUNT )
    GRAVITY IS FREE FROM ABERRATION  (-0.0057 DEG.)

    TIME : (JED-2451545.0)/36525     EPOCH = J2000.0
    """
    B0=36000.7695
    C0=280.4659
    A=np.array([19147e-4, 200e-4, 48e-4, 20e-4, 18e-4, 18e-4,\
        15e-4,  13e-4,  7e-4,  7e-4,  7e-4,  6e-4,\
        5e-4,   5e-4,  4e-4,  4e-4])
    B=np.array([35999.050, 71998.1,  1934, 32964,    19,\
        445267   , 45038 , 22519, 65929,  3035,\
        9038   , 33718 ,   155,  2281, 29930,\
        31557   ])
    C=np.array([267.520, 265.1, 145, 158, 159, 208,\
        254.   , 352 ,  45, 110,  64, 316,\
        118.   , 221 ,  48, 161])
    RAD=0.0174532925199433
    A[0] = 1.9147-0.0048*time
    tempb=(B*time+C)*RAD
    amp=A*np.cos(tempb)
    sunlon=np.sum(amp)
    sunlon=(sunlon+B0*time+C0)*RAD
    return sunlon    
    
def moonlongitude(time):
    """
    APPARENT LONGITUDE OF THE MOON IN RADIANS
    
    TIME : (JEe-2451545.0)/36525          EPOCH = J2000.0
    """    
    B0=481267.8809
    C0=218.3162
    A=np.array([62888.e-4, 12740.e-4, 6583.e-4, 2136.e-4, 1851.e-4,\
        1144.e-4, 588.e-4, 571.e-4, 533.e-4, 458.e-4, 409.e-4,\
        347.e-4, 304.e-4, 154.e-4, 125.e-4, 110.e-4, 107.e-4,\
        100.e-4, 85.e-4, 79.e-4, 68.e-4, 52.e-4, 50.e-4, 40.e-4,\
        40.e-4, 40.e-4, 38.e-4, 37.e-4, 28.e-4, 27.e-4, 26.e-4,\
        24.e-4, 23.e-4, 22.e-4, 21.e-4, 21.e-4, 21.e-4, 18.e-4,\
        16.e-4, 12.e-4, 11.e-4,  9.e-4,  8.e-4,  7.e-4,  7.e-4,\
        7.e-4,  7.e-4,  6.e-4,  6.e-4,  5.e-4,  5.e-4,  5.e-4,\
        4.e-4,  4.e-4,  3.e-4,  3.e-4,  3.e-4,  3.e-4,  3.e-4,\
        3.e-4,  3.e-4])
    B=np.array([477198.868,  413335.35, 890534.22, 954397.74,\
        35999.05 ,  966404.0 ,  63863.5 , 377336.3 ,\
        1367733.1  ,  854535.2 , 441199.8 , 445267.1 ,\
        513197.9, 75870,1443603, 489205,1303870,\
        1431597, 826671, 449334, 926533,  31932,\
        481266,1331734,1844932,    133,1781068,\
        541062,   1934, 918399,1379739,  99863,\
        922466, 818536, 990397,  71998, 341337,\
        401329,1856938,1267871,1920802, 858602,\
        1403732, 790672, 405201, 485333,  27864,\
        111869,2258267,1908795,1745069, 509131,\
        39871,  12006, 958465, 381404, 349472,\
        1808933, 549197,   4067,2322131.])
    C=np.array([44.963, 10.74, 145.70, 179.93, 87.53,  276.5,\
        124.2, 13.2, 280.7, 148.2, 47.4, 27.9, 222.5,\
        41,  52, 142, 246, 315, 111, 188,\
        323, 107, 205, 283,  56,  29,  21,\
        259, 145, 182,  17, 122, 163, 151,\
        357,  85,  16, 274, 152, 249, 186,\
        129,  98, 114,  50, 186, 127,  38,\
        156,  90,  24, 242, 223, 187, 340,\
        354, 337,  58, 220,  70, 191])
    RAD=0.0174532925199433
    tempb=(B*time+C)*RAD
    amp=A*np.cos(tempb)
    moonlon=np.sum(amp)
    moonlon=(moonlon+B0*time+C0)*RAD 
    return moonlon
    
def gravityeffect(h0,s0,p0,tl,amp,phases):
    """
    purpose:  To compute the gravity effect due to each        
    tidal component.
    """
    dtor = np.pi/180
    R90 = np.pi/2       
    R360 = np.pi*2
    p0 = np.remainder(p0,R360)
    h0 = np.remainder(h0,R360)
    s0 = np.remainder(s0,R360)
    tl = np.remainder(tl,R360)
    arg=np.zeros(11)
    # argument                     component
    arg[0] = 2*tl - 2*s0          # M2
    arg[1] = 2*tl - 2*h0          # S2
    arg[2] = tl - R90             # K1
    arg[3] = tl - 2*s0 + R90      # O1
    arg[4] = 2*tl - 3*s0 + p0     # N2
    arg[5] = tl - 2*h0 + R90      # P1
    arg[6] = 2*tl                 # K2
    arg[7] = tl - 3*s0 + p0 +R90  # Q1
    arg[8] = 2*s0                 # Mf
    arg[9] = s0 - p0              # Mm
    arg[10] = 2*h0                # Ssa
    totaleffect=np.sum([amp[i]*np.cos(arg[i] - phases[i]*dtor) for i in range(len(arg))])
    return totaleffect    