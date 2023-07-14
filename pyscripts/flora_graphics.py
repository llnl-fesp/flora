
import matplotlib.pyplot as plt
from flora import *
import numpy as np


def plots1():
    iz1 = 1
    iz2 = np.argmin(np.abs(glr.z - 140*glr.zmax/300))
    iz3 = np.argmin(np.abs(glr.z - 150*glr.zmax/300))
    iz4 = np.argmin(np.abs(glr.z - 270*glr.zmax/300))
    iz5 = np.argmax(glr.z)

    xl = 'Z [cm]'
    
    fig,axs = plt.subplots(2)
    fig.suptitle(utl.probid.tostring().strip().decode())
    axs[0].plot(glr.z[iz2:iz5],glr.pperp[iz2:iz5,1])
    axs[0].plot(glr.z[iz2:iz5],glr.pperpold[iz2:iz5],linestyle='dashed')
    axs[0].set(xlabel=xl,ylabel='pperp,pperpold')
    axs[0].label_outer()
    axs[1].plot(glr.z[iz2:iz5],glr.pperp[iz2:iz5,1])
    axs[1].plot(glr.z[iz2:iz5],glr.ppar[iz2:iz5,1])
    axs[1].set(xlabel=xl,ylabel='pperp,ppar')
    axs[1].label_outer()
    plt.show()

    fig,axs = plt.subplots(2)
    fig.suptitle(utl.probid.tostring().strip().decode())
    axs[0].plot(glr.z[iz2:iz5],glr.pperpold[iz2:iz5],linestyle='dashed')
    axs[0].plot(glr.z[iz2:iz5],glr.pparold[iz2:iz5],linestyle='dashed')
    axs[0].set(xlabel=xl,ylabel='pperpold,pparold')
    axs[0].label_outer()
    axs[1].plot(glr.z[iz2:iz5],glr.ppar[iz2:iz5,1])
    axs[1].plot(glr.z[iz2:iz5],glr.pparold[iz2:iz5],linestyle='dashed')
    axs[1].set(xlabel=xl,ylabel='ppar,pparold')
    axs[1].label_outer()
    plt.show()

    fig,axs = plt.subplots(2)
    fig.suptitle(utl.probid.tostring().strip().decode())
    axs[0].plot(glr.bvac[iz3:iz5]/glr.bvac[139],glr.pperp[iz3:iz5,1],linestyle='dashed')
    axs[0].plot(glr.bvac[iz3:iz5]/glr.bvac[139],glr.pperpold[iz3:iz5],linestyle='dashed')
    axs[0].set(xlabel='bvac/bvacinj',ylabel='pperp,pperpold')
    axs[1].axis('off')
    plt.show()
