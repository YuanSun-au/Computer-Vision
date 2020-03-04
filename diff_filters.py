import numpy as np
from scipy import signal, misc, ndimage, stats
from scipy.ndimage.filters import convolve
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
# from mpl_toolkits.mplot3d import Axes3D

def gaussian_filter( img, sigma=1, precision=3, h=1):

    # calculate gaussian kernel
    length = np.max( [1, np.ceil(precision*sigma/h).astype('int') ])
    gaussian1D = signal.gaussian(M=length*2+1, std=sigma/h)
    gaussian2D = np.outer(gaussian1D, gaussian1D) / (2*np.pi*(sigma**2))
    gaussian2D /= gaussian2D.sum()

    # # calculate gaussian kernel
    # length = np.max( [1, np.ceil(precision*sigma/h).astype('int') ])
    # wid = np.arange(-length, length+1)*h
    # x = np.array([ [ [j,i] for j in wid] for i in wid])
    # kernel = stats.multivariate_normal.pdf(x, mean=2*[0], cov=sigma**2*np.eye(2))
    # kernel /= kernel.sum()

    # fig, (ax1,ax2) = plt.subplots(1,2)
    # im = ax1.imshow(kernel)
    # cbar = ax1.figure.colorbar(im, ax=ax1)
    # im = ax2.imshow(gaussian2D)
    # cbar = ax2.figure.colorbar(im, ax=ax2)
    # plt.show()

    # return ndimage.gaussian_filter(img, sigma=sigma, mode='constant')
    return signal.convolve2d(img, gaussian2D, boundary='fill', mode='same')
    
def EA(img, D, tau=0.1, hx=1, hy=1):
    
    imgp = np.pad(img,((1,1),(1,1)), mode='symmetric') # reflect
    Dp = np.pad(D, ((1,1),(1,1),(0,0),(0,0)), mode='constant')
    
    ny,nx = img.shape
    
    stencil = lambda D,j,i: \
        np.array([ [
                        (-D[j,i-1,0,1]-D[j+1,i,0,1])/(4*hx*hy),
                        (+D[j+1,i,1,1]+D[j,i,1,1])/(2*hy*hy),
                        (+D[j,i+1,0,1]+D[j+1,i,0,1])/(4*hx*hy)
                    ] ,
                    [
                        (+D[j,i-1,0,0]+D[j,i,0,0])/(2*hx*hx),
                       -(+D[j,i-1,0,0]+2*D[j,i,0,0]+D[j,i+1,0,0])/(2*hx*hx) -(+D[j-1,i,1,1]+2*D[j,i,1,1]+D[j+1,i,1,1])/(2*hy*hy),
                        (+D[j,i+1,0,0]+D[j,i,0,0])/(2*hx*hx)
                    ],
                    [
                        (+D[j,i-1,0,1]+D[j-1,i,0,1])/(4*hx*hy),
                        (+D[j-1,i,1,1]+D[j,i,1,1])/(2*hy*hy),
                        (-D[j,i+1,0,1]-D[j-1,i,0,1])/(4*hx*hy)
                    ]
                ])
    imgp_new = np.copy(imgp)
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            imgp_new[j,i] = imgp[j,i] + tau * np.sum( np.flipud(stencil(Dp,j,i)) * imgp[j-1:j+2,i-1:i+2] )
    return imgp_new[1:-1,1:-1]


def CEAD(img, sigma=1, tau=0.1, tf=1, hx=1, hy=1):
    ''' Coherence-Enhancing Anysotropic Diffusion'''
    pass
    # return signal.convolve2d(img, kernel, boundary='fill', mode='same')


def EEAD(u, sigma=1, lmb=3, tau=0.1, tf=1, hx=1, hy=1):
    ''' Edge-Enhancing Anysotropic Diffusion'''

    ny,nx = u.shape

    dx = np.array([[1,0,-1]])/(2*hx)
    dy = np.array([[1,0,-1]]).T/(2*hy)

    g = lambda s2: 1/(1+s2/lmb**2)

    for t in np.arange(0,tf,tau):
        print(t)

        # calculate gradient tensor
        u_s = ndimage.gaussian_filter(u, sigma=sigma, mode='constant')
        du = np.zeros([ny,nx,2])
        du[:,:,0] = convolve(u_s, dx, mode='mirror')
        du[:,:,1] = convolve(u_s, dy, mode='mirror')

        # plt.imshow( du[:,:,0] ); plt.show()
        # plt.imshow( du[:,:,1] ); plt.show()
        # plt.imshow( u ); plt.show()

        # calculate structure tensor
        # J = np.array([ [ np.outer(du[j,i,:],du[j,i,:]) for j in range(n[1])] for i in range(n[0])])
        # for i in range(2):
        #     for j in range(2):
        #         J[:,:,j,i] = ndimage.gaussian_filter(J[:,:,j,i], sigma=sigma, mode='constant')
        
        # calculate difusion tensor
        D = np.zeros([ny,nx,2,2])
        for i in range(nx):
            for j in range(ny):
                # eigval, eigvec = np.linalg.eig( J[j,i,:,:] )
                eigval = [g(np.linalg.norm(du[j,i,:])**2),1]
                eigvec = np.zeros([2,2])
                eigvec[:,0] = du[j,i,:]
                if np.linalg.norm(eigvec[:,0]) > 1e-10:
                    eigvec[:,0] /= np.linalg.norm(eigvec[:,0])
                else:
                    eigvec[:,0] = 0
                eigvec[:,1] = [[0,1],[-1,0]] @ eigvec[:,0]
                D[j,i,:,:] = eigvec @ np.diag(eigval) @ eigvec.T
                # D[j,i,:,:] = np.eye(2)

        # perform one diffusion step
        u = EA(u,D,tau,hx,hy)

    return u


''' MAIN '''
sigma = 0.5
precision = 3
h = 1

img = mpimg.imread('./cv19/images/sbrain.pgm')[::2,::2]
# img = misc.face(gray=True).astype(np.float32)
# img_g1 = gaussian_filter(img, sigma=sigma, precision=precision, h=h)
img_g1 = ndimage.gaussian_filter(img, sigma=sigma, mode='constant')
img_g2 = EEAD(img, sigma=1, tau=0.2, tf=10, hx=h, hy=h)
img_g3 = EEAD(img_g2, sigma=1, tau=0.2, tf=10, hx=h, hy=h)

fig, ax = plt.subplots(1,4, figsize=(12, 4)) # , figsize=(2,2)
ax[0].imshow(img, cmap='gray')
ax[0].set_title('Original')
ax[0].set_axis_off()
ax[1].imshow(img_g1, cmap='gray')
ax[1].set_title('Blured image my function')
ax[1].set_axis_off()
ax[2].imshow(img_g2, cmap='gray')
ax[2].set_title('Blured image standard lib')
ax[2].set_axis_off()
ax[3].imshow(img_g3, cmap='gray')
ax[3].set_title('Blured image standard lib')
ax[3].set_axis_off()
plt.show()
