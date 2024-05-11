"""
tests.py
Chase Funkhouser

This script supplies three different code tests as follows:

    Test 1: Use an input namelist with all force parameters set to zero. The code should then find the critical radius at the parker radius.

    Test 2: Use an input namelist with a constant force term. This can be solved analytically.

    Test 3: Use an input namelist with a full polynomial force term. Ensure that the mass loss rate is constant over the full domain.

This script can either be called through the shell script test.sh, or can be called using the command line interface.
To run, use the -t flag followed by the test number to run.
Be sure that the code was previously run using the correct test input namelist.

"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

KMS_TO_CMS=1.e5
EPS=1.e-4
A0=1.e-5


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


def test1(cs,gm,rcrit,eps):
    """ Determines the truth condition of the first test.
    The first test compares the output from a model with no force to the original Parker solution."""

    r_parker=gm/2/(cs*KMS_TO_CMS)**2
    print('Parker radius=',r_parker)
    print('Critical radius=',rcrit)
    error = abs(rcrit-r_parker)/r_parker
    print('Relative error=',error)
    if error < eps:
        return True
    return False


def test2(cs,gm,rcrit,g0,eps):
    """ Returns the truth condition of the second test.
    This test compares the output from a constant acceleration model to the derived solution."""

    cs2 = (cs*KMS_TO_CMS)**2
    rtrue = (-2*cs2+np.sqrt(4*cs2**2+4*g0*gm))/(2*g0)
    print('Analytic critical radius=',rtrue)
    print('Critical radius=',rcrit)
    error = abs(rcrit-rtrue)/rtrue
    print('Relative error=',error)
    if error < eps:
        return True
    return False


def test3(mdot,eps):
    """ Returns the truth condition of the third test.
    This test ensures that the mdot value is constant over radius, even for large forces."""

    maxmdot = np.max(mdot)
    minmdot = np.min(mdot)
    error = abs(maxmdot-minmdot)/minmdot
    if error < eps:
        return True
    return False


def main():
    parser=argparse.ArgumentParser()

    #-t TEST_TYPE
    parser.add_argument('-t',dest='test',help='Test type')
    args = parser.parse_args()
    test = args.test
    print('Test number: '+test)

    nrad,temp,adia,mplan,cs,gm,rcrit = read_setup_file('setup.data')

    match (test):
        case '1':
            nrad,rad,force,vel=read_vel_file('final_velocity.data')
            t1=test1(cs,gm,rcrit,EPS)
            if t1:
                print('Test 1 successful within error.')
            else:
                print('Test 1 failed.')

        case '2':
            nrad,rad,force,vel=read_vel_file('final_velocity.data')
            t2=test2(cs,gm,rcrit,A0,EPS)
            if t2:
                print('Test 2 successful within error.')
            else:
                print('Test 2 failed.')

        case '3':
            nrad,rad,density,mdot=read_density_file('density_mdot.data')
            t3=test3(mdot,EPS)
            if t3:
                print('Test 3 successful within error.')
            else:
                print('Test 3 failed.')

        case _:
            raise ValueError("Not a supported plot type.")


if __name__ == '__main__':
    main()
