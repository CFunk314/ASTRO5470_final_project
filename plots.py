"""
plots.py
Chase Funkhouser
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

KMS_TO_CMS=1.e5


def read_setup_file(filename):
    """ Read setup.data file for a run and return parameters. """

    data=np.loadtxt(filename,usecols=1)
    
    nrad=int(data[0])
    temp=data[3]
    adia=data[4]
    mplan=data[5]
    cs=data[6]
    gm=data[7]

    return nrad,temp,adia,mplan,cs,gm


def read_vel_file(filename):
    """ Read radius and velocity data from a .data file as arrays. """

    f=open(filename,'r')
    h=f.readline()
    h=h.split()
    nrad=int(h[1])
    f.close()

    data=np.genfromtxt(filename,names=True,skip_header=1)
    rad=data['rad']
    force=data['force']
    vel=data['vel']

    assert len(rad) == nrad
    assert len(force) == nrad
    assert len(vel) == nrad

    return nrad,rad,force,vel


def read_density_file(filename):
    """ Read density and mdot data from a .data file as arrays. """

    f=open(filename,'r')
    h=f.readline()
    h=h.split()
    nrad=int(h[1])
    f.close()

    data=np.genfromtxt(filename,names=True,skip_header=1)
    rad=data['rad']
    density=data['density']
    mdot=data['mdot']

    assert len(rad) == nrad
    assert len(density) == nrad
    assert len(mdot) == nrad

    return nrad,rad,density,mdot


def plot_vel(nrad,rad,vel,cs,gm,filename):
    """ Plot velocity throughout atmosphere. """

    plt.plot(rad,vel,'b-',label='velocity')

    plt.axhline(y=cs,c='r',linestyle='-',label='sound speed')

    parker=gm/2/(cs*KMS_TO_CMS)**2
    plt.axvline(x=parker,c='b',linestyle='dotted',label='parker point')

    plt.ylim(-0.5,20)
    plt.xscale('log')
    plt.legend(loc='best')
    plt.xlabel(r'$r\ [cm]$')
    plt.ylabel(r'$v\ [km/s]$')
    plt.title('velocity, nrad='+str(nrad))
    plt.savefig(filename,format='png')
    plt.show()
    plt.close()


def plot_dens(nrad,rad,density,filename):
    """ Plot density throughout atmosphere. """

    plt.plot(rad,density,'b-',label='density')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.xlabel(r'$r\ [cm]$')
    plt.ylabel(r'$\rho\ [g/cm^3]$')
    plt.title('density, nrad='+str(nrad))
    plt.savefig(filename,format='png')
    plt.show()
    plt.close()


def plot_mdot(nrad,rad,mdot,filename):
    """ Plot mdot throughout atmosphere. """

    plt.plot(rad,mdot,'b-',label='mdot')

    plt.xscale('log')
    plt.legend(loc='best')
    plt.xlabel(r'$r\ [cm]$')
    plt.ylabel(r'$\dot{M}\ [g/s]$')
    plt.title('mdot, nrad='+str(nrad))
    plt.savefig(filename,format='png')
    plt.show()
    plt.close()


def main():
    parser=argparse.ArgumentParser()

    #-i INFILE, -p PLOT_TYPE
    parser.add_argument('-i',dest='infile',help='Infile')
    parser.add_argument('-p',dest='plot',help='Plot type')
    args = parser.parse_args()
    infile = args.infile
    plot_type=args.plot
    print('Infile: '+infile)
    print('Plotting: '+plot_type)
    outfile = infile[:-4]+'png' 
    print('Outfile: '+outfile)

    nrad,temp,adia,mplan,cs,gm = read_setup_file('setup.data')

    match (plot_type):
        case 'velocity':
            nrad,rad,force,vel=read_vel_file(infile)
            plot_vel(nrad,rad,vel,cs,gm,'velocity_plot.png')
        case 'density':
            nrad,rad,density,mdot=read_density_file(infile)
            plot_dens(nrad,rad,density,'density_plot.png')
        case 'mdot':
            nrad,rad,density,mdot=read_density_file(infile)
            plot_mdot(nrad,rad,mdot,'mdot_plot.png')
        case _:
            raise ValueError("Not a supported plot type")


if __name__ == '__main__':
    main()
