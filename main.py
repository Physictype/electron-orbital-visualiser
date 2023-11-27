import numpy as np
import os
from tqdm import tqdm
import math
import time

def associated_laguerre(p,q):
    return ((-1)**p)*laguerre(p+q).deriv(p)
def laguerre(q):
    res = []
    for k in range(q,-1,-1):
        res.append(math.comb(q,k)*(-1)**k/math.factorial(k))
    return np.poly1d(res)


def associated_legendre(m,l):
    return lambda x: (-1)**m*np.sqrt(squared_associated_legendre(m,l)(x))
def squared_associated_legendre(m,l):
    return np.poly1d([-1,0,1])**m*(legendre(l).deriv(m))**2
def legendre(l):
    return (1/(2**l*math.factorial(l)))*(np.poly1d([1,0,-1])**l).deriv(l)

def spherical_harmonic(m,l):
    return lambda theta,phi: np.sqrt((2*l+1)/(4*np.pi)*math.factorial(l-m)/math.factorial(l+m))*np.exp(1j*m*phi)*associated_legendre(m,l)(np.cos(theta))

def r(x,y,z):
    return np.sqrt(x**2+y**2+z**2)
def theta(x,y,z):
    return np.arctan(y/x)
def phi(x,y,z):
    return np.arccos(z/r(x,y,z))
def q(x,y,z):
    return Z*r(x,y,z)/a0/n

Z = 1
a0 = 1

def wave_func(n,l,m):
    if l < 0 or l != np.round(l) or l >= n:
        raise ValueError(f"l={l} must be a nonnegative integer from 0 to {n-1}=n-1.")
    if np.abs(m) > l or m != np.round(m):
        raise ValueError(f"m={m} must be an integer between -l={-l} and {l}=l.")
    try:
        return lambda x,y,z: np.sqrt(np.power(2/n/a0,3)*math.factorial(n-l-1)/(2*n)/math.factorial(n+l))*np.exp(-q(x,y,z))*np.power(2*q(x,y,z),l)*associated_laguerre(2*l+1,n-l-1)(2*q(x,y,z))*spherical_harmonic(m,l)(theta(x,y,z),phi(x,y,z))
    except ZeroDivisionError:
        return 0

os.system("clear")
input("Welcome to the electron orbital renderer.\nPress [Enter] to continue.")
os.system("clear")
input("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\nPress [Enter] to continue.")
os.system("clear")
print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")

print("BOUNDS")
u_min = float(input("Enter the minimum u value: ") or -12)
u_max = float(input("Enter the maximum u value: ") or 12)
v_min = float(input("Enter the minimum v value: ") or -12)
v_max = float(input("Enter the maximum v value: ") or 12)
w = float(input("Enter the value of w for where to slice: ") or 0)
print("\n")
wave_min = float(input("Enter the minimum wave function value: ") or 0)
wave_max = float(input("Enter the maximum wave function value: ") or 0)

os.system("clear")
print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")
print("PLANE ANGLES")
a_x = input("Would you like to have variable x rotation? (y if yes, anything else for no) ") == "y"
if a_x:
    x_rot_start = float(input("\tWhat would you like to have the x rotation start at? (in degrees for simplicity) ") or 0) * np.pi/180
    x_rot_end = float(input("\tWhat would you like to have the x rotation end at? ") or 360) * np.pi/180
else:
    x_rot = float(input("\tHow much would you like the plane rotated around the x axis? (in degrees for simplicity) ") or 0) * np.pi/180
    x_rot_start = x_rot
    x_rot_end = x_rot
a_y = input("Would you like to have variable y rotation? ") == "y"
if a_y:
    y_rot_start = float(input("\tWhat would you like to have the y rotation start at? (in degrees for simplicity) ") or 0) * np.pi/180
    y_rot_end = float(input("\tWhat would you like to have the y rotation end at? ") or 360) * np.pi/180
else:
    y_rot = float(input("\tHow much would you like the plane rotated around the y axis? (in degrees for simplicity) ") or 0) * np.pi/180
    y_rot_start = y_rot
    y_rot_end = y_rot
a_z = input("Would you like to have variable z rotation? ") == "y"
if a_z:
    z_rot_start = float(input("\tWhat would you like to have the z rotation start at? (in degrees for simplicity) ") or 0) * np.pi/180
    z_rot_end = float(input("\tWhat would you like to have the z rotation end at? ") or 360) * np.pi/180
else:
    z_rot = float(input("\tHow much would you like the plane rotated around the z axis? (in degrees for simplicity) ") or 0) * np.pi/180
    z_rot_start = z_rot
    z_rot_end = z_rot

if a_x or a_y or a_z:
    rot_time = int(input("How long would you like the animation to last in frames? ") or 60)
    
# plane_theta
# plane_phi

while True:
    os.system("clear")
    print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")
    print("CONFIGURATION*")
    try:
        n = int(input("Enter the electron's n: "))
    except:
        while True:
            input("Please enter a number.")
            os.system("clear")
            print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")
            print("CONFIGURATION*")
            try:
                n = int(input("Enter the electron's n: "))
                break
            except:
                continue
    try:
        l = int(input("Enter the electron's l: "))
    except:
        while True:
            input("Please enter a number.")
            os.system("clear")
            print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")
            print("CONFIGURATION*")
            print(f"Enter the electron's n: {n}")
            try:
                l = int(input("Enter the electron's l: "))
                break
            except:
                continue
    try:
        m = int(input("Enter the electron's m: "))
    except:
        while True:
            input("Please enter a number.")
            os.system("clear")
            print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")
            print("CONFIGURATION*")
            print(f"Enter the electron's n: {n}")
            print(f"Enter the electron's l: {l}")
            try:
                m = int(input("Enter the electron's m: "))
                break
            except:
                continue
    try:
        wave_func(n,l,m)
        break
    except ValueError as e:
        input(e)

os.system("clear")
print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")
print("IMAGE DIMENSIONS")
width = int(input("Enter the width in pixels: ") or 256)
height = int(input("Enter the height in pixels: ") or 256)

os.system("clear")
print("Welcome to the electron orbital renderer.\nBelow there will be configuration settings. Only settings or categories marked with a * are needed, others can be defaulted.\n\n")
print("MISCELLANEOUS")
filepath = input("Enter file path: ") or "output.png"

os.system("clear")
print("Loading...")

wf = wave_func(n,l,m)

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)

    b = -axis[0] * np.sin(theta / 2.0)
    c = -axis[1] * np.sin(theta / 2.0)
    d = -axis[2] * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc - ad), 2 * (bd + ac)],
                     [2 * (bc + ad), aa + cc - bb - dd, 2 * (cd - ab)],
                     [2 * (bd - ac), 2 * (cd + ab), aa + dd - bb - cc]])

def vrotate(v,axis,theta):
    return np.dot(v,rotation_matrix(axis, theta))

i = [1,0,0]
j = [0,1,0]
k = [0,0,1]
from PIL import Image
def generate_img(ri,rj,rk):
    data = np.zeros((width, height, 3))
    u,v = np.meshgrid(np.arange(width),np.arange(height))
    u = u.flatten()
    v = v.flatten()
    scaled_u = (u/(width-1))*(u_max-u_min)+u_min
    scaled_v = (v/(width-1))*(v_max-v_min)+v_min
    vec = ri*scaled_u[:,np.newaxis]+rj*scaled_v[:,np.newaxis]+rk*w
    vec = np.transpose(vec)
    if wave_max == 0:
        c = np.power(np.abs(wf(vec[0],vec[1],vec[2])),2)*1000
    else:
        c = (np.abs(wf(vec[0],vec[1],vec[2]))**2-wave_min)/(wave_max-wave_min)*255
    c = c.reshape((width,height))
    data=np.array([c,c,c]).swapaxes(0,2).swapaxes(0,1)
    if wave_max == 0:
        data = (data-wave_min)/(np.max(data)-wave_min)*255
    data = data.astype(np.uint8)
    return Image.fromarray(data)

if not (a_x or a_y or a_z):
    rot_time = 0
imgs = []
for t in tqdm(range(rot_time)):
    x_rot = x_rot_start+t/(rot_time-1)*(x_rot_end-x_rot_start)
    y_rot = y_rot_start+t/(rot_time-1)*(y_rot_end-y_rot_start)
    z_rot = z_rot_start+t/(rot_time-1)*(z_rot_end-z_rot_start)
    rotated_i = vrotate(vrotate(vrotate(i,[1,0,0],x_rot),[0,1,0],y_rot),[0,0,1],z_rot)
    rotated_j = vrotate(vrotate(vrotate(j,[1,0,0],x_rot),[0,1,0],y_rot),[0,0,1],z_rot)
    rotated_k = vrotate(vrotate(vrotate(k,[1,0,0],x_rot),[0,1,0],y_rot),[0,0,1],z_rot)
    imgs.append(generate_img(rotated_i,rotated_j,rotated_k))
imgs[0].save(filepath,save_all=True,append_images=imgs[1:])
print(f"Done. Check the image in {filepath}")