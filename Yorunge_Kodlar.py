"""
Kadir Can Kızılarslan
20183405033
Target Body: 33 Polyhymnia 
"""

import numpy as np
import math
def derece(x):
 x = 180/math.pi
 return x
X = -1.928205198918929E+00
Y =-2.306271811098207E+00
Z = -6.478259352027139E-02
VX = 9.067442195237805E-03
VY =-3.388319970653067E-03 
VZ = -1.546560685584234E-04
k = 2.9591220828411951E-04
def r_(r):
 r = np.array([X,Y,Z])
 return r
def v_(v):
 v = np.array([VX,VY,VZ])
 return v
def v(v):
 v = np.sqrt(VX**2+ VY**2 + VZ**2)
 return v
r_vector= np.array([X,Y,Z])
v_vector=np.array([VX,VY,VZ])
r = np.sqrt((X**2 + Y**2 + Z**2))
print("r= ",r)
v = v("v")
print("v= ",v)
print("r.v ve r x v= ",np.dot(r_("r"),v_("v")), np.cross(r_("r"),v_("v")))
a = 2/r - (v**2/k)
print("bulduğum a= ",1/a)
a = round(1/a,15)
a_gercek = 2.870059095472908E+00
print("a mutlak hata =", np.abs(a_gercek - a))
# Hata nedeni bulduğum değerin son basamağı horizonda verilen değerden bir basmak fazla olarak çıkmasından kaynaklıydı o yüzden 15 basamağa yuvarlayarak işlemi yaptım.
# Yuvarlamadığım takdirde aldığım değer = 2.8700590954729086
def ecc(r_vector,v_vector,k):
    return (((v**2/k)-(1.0/r))*r_("r")-(np.dot(r_("r"),v_("v"))/(k))*v_("v"))
def ecc_1(r_vec,v_vec,k):
    return (np.cross(v_("v"),np.cross(r_("r"),v_("v")))/k - r_("r")/r)
ecc_vec = ecc(r_("r"),v_("v"),k)
eccentricity=np.linalg.norm(ecc_vec)
print("Basıklık vektörü =",ecc_vec,ecc_1(r_("r"),v_("v"),k))
print("Basıklık =",eccentricity,np.linalg.norm(ecc_1(r_vector,v_vector,k)))
EC = 3.348653410238699E-01
EC_error = EC - 0.3348653410244694
print("Hata=", np.abs(EC_error))
h_vec= np.cross(r_vector,v_vector)
print("h_vec=",h_vec)
h= np.sqrt(0.00013717**2 + (-0.00088562**2) + 0.02744536**2)
def h_(h):
 h = np.sqrt(0.00013717**2 + (-0.00088562**2) + 0.02744536**2)
 return h
print("h = ",h)
i_deg=derece("x")*(np.arccos(h_vec[2]/np.linalg.norm(h_vec)))
IN = 1.870231780893287E+00
print("i=",i_deg,"err=",IN-i_deg)
N_vec=np.array([-h_vec[1],h_vec[0],0])
print("N_vec:",N_vec)
#Büyük Omega
Omega= derece("x")*(np.arccos((-h_vec[1]/np.linalg.norm(N_vec))))
OM= 8.804597359125099E+00
print("Omega=",Omega,"err=",OM-Omega)
#Küçük omega
omega = 360-(derece("x")*(np.arccos(np.dot(N_vec,ecc_vec)/(np.linalg.norm(N_vec)*np.linalg.norm(ecc_vec)))))
W = 3.384082266251065E+02
print("omega=",omega,"err=",W-omega)
#Öz kiriş
rc= (np.linalg.norm(h_vec)**2)/k
# rc = a*(1-eccentricity**2)
print("rc err kontrol=",rc -(a*(1 - eccentricity**2)))
#x
x = (rc-r)/eccentricity
print("x=",x,"AB")
y = np.dot(r_vector,v_vector)*np.sqrt(rc/k)/eccentricity
print("y=",y,"AB")
print("r err=",np.sqrt(x**2+y**2)-r)
cos_E = x/a+eccentricity
sin_E = y/(a*np.sqrt(1-eccentricity**2))
print("cos(E)=",cos_E,"rad sin(E)=",sin_E,"rad")
if np.abs(cos_E) <= 0.707107 :
 E = np.arccos(np.abs(cos_E))
if cos_E<0 and sin_E < 0:
 E = np.pi+E
print("E=",derece("x")*(E))
M= E-eccentricity*np.sin(E)
print("M=",derece("x")*(M),"derece")
MA = 2.808092426023672E+02
print("M err=",MA-derece("x")*(M))
