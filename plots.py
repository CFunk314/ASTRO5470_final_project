"""
plots.py
Chase Funkhouser

This script supplies various plotting functionality including:

    velocity plot

    density plot

    mdot plot

This script is called with two command-line arguments:

    -i specifies the input file

    -p specifies the plot type (velocity, density, mdot)

Example of calling the script:

    python plots.py -i final_velocity.data -p velocity

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
    rcrit=data[8]

    return nrad,temp,adia,mplan,cs,gm,rcrit


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


def plot_vel(nrad,rad,vel,cs,gm,rcrit,filename):
    """ Plot velocity throughout atmosphere. """
    
    # plot velocity
    plt.plot(rad,vel,'b-',label='velocity')

    # plot the sound speed
    plt.axhline(y=cs,c='r',linestyle='-',label='sound speed')

    # calculate and plot the Parker radius given sound speed and GM
    parker=gm/2/(cs*KMS_TO_CMS)**2
    plt.axvline(x=parker,c='g',linestyle='dashed',label='parker point')

    # plot the actual critical radius
    plt.axvline(x=rcrit,c='b',linestyle='dotted',label='critical point')

    plt.ylim(-0.5,20)
    plt.xscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel(r'$r\ [cm]$')
    plt.ylabel(r'$v\ [km/s]$')
    plt.title('velocity, nrad='+str(nrad))
    plt.savefig(filename,format='png')
    plt.show()
    plt.close()


def plot_dens(nrad,rad,density,rcrit,filename):
    """ Plot density throughout atmosphere. """

    # calculate the closest index to the critical radius
    diff=np.inf
    indx=0
    for i in range(nrad):
        newdiff=abs(rcrit-rad[i])
        if newdiff < diff:
            diff=newdiff
            indx=i
    
    # plot density
    plt.plot(rad,density,'b-',label='density')
    
    # plot critical radius
    plt.axvline(x=rcrit,c='b',linestyle='dotted',label='critical point')

    # plot density at the critical radius (critical density)
    plt.axhline(y=density[indx],c='r',linestyle='-',label='critical density')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel(r'$r\ [cm]$')
    plt.ylabel(r'$\rho\ [g/cm^3]$')
    plt.title('density, nrad='+str(nrad))
    plt.savefig(filename,format='png')
    plt.show()
    plt.close()


def plot_mdot(nrad,rad,mdot,filename):
    """ Plot mdot throughout atmosphere. """

    # plot mdot
    plt.plot(rad,mdot,'b-',label='mdot')

    plt.xscale('log')
    plt.grid()
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

    nrad,temp,adia,mplan,cs,gm,rcrit = read_setup_file('setup.data')

    match (plot_type):
        case 'velocity':
            nrad,rad,force,vel=read_vel_file(infile)
            plot_vel(nrad,rad,vel,cs,gm,rcrit,'velocity_plot.png')
        case 'density':
            nrad,rad,density,mdot=read_density_file(infile)
            plot_dens(nrad,rad,density,rcrit,'density_plot.png')
        case 'mdot':
            nrad,rad,density,mdot=read_density_file(infile)
            plot_mdot(nrad,rad,mdot,'mdot_plot.png')
        case _:
            raise ValueError("Not a supported plot type")


if __name__ == '__main__':
    main()
