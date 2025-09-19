# 558 Question 2, Andrew Mckenna 

import numpy as np 
from matplotlib import pyplot as plt

try:
    from matplotlib import rc
    rc( 'text', usetex = True ) # use nice LaTeX font
except:
    print('error with latex font')


plt.ion() # interactive plotting 


V0 = 1.0 # depth of the well in Ry
xres = 400 
Rres = 50
initial_well_buffer = 0.2 # how far appart are the well in units of a0


# define functions:
def sorted_eig( Hamiltonian ): # by default eig does not sort eigenvalues:

    eigenvalues, eigenvectors = np.linalg.eig( Hamiltonian ) 
    idx = np.argsort( eigenvalues )
    eps_sorted = eigenvalues[idx]
    psi_sorted = eigenvectors[:, idx]

    return eps_sorted, psi_sorted

def hc( matrix ): # computes the Hermitian conjugate
    return np.transpose(np.conj( matrix ))


xvals = np.linspace( -12.0, 12.0, xres )
dx = xvals[1] - xvals[0]

well_inc = np.ceil(0.5*(1 + initial_well_buffer)/dx).astype(int)
Rvals = np.zeros((Rres, ))

Evals = np.zeros((Rres, 2))
Evals_TBA = np.zeros((Rres, 2))

# define the kinetic part of the Hamiltonian
KE = ( -np.diag( np.ones((xres - 1, )), -1 ) 
    + 2*np.eye( xres ) 
    - np.diag( np.ones((xres - 1, )), 1) )/(dx**2);

# get the potential of a single well 
Vsingle_well = np.zeros((xres, ))
for i in range(xres): # check if we are in the well:
    if -0.5 < xvals[i] and xvals[i] < 0.5:
        Vsingle_well[i] = -1.0


# initialize the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize = (8, 6))


for i in range(Rres):

    Vdouble_well = np.zeros((xres, ))
    
    Vleft_well = np.roll(Vsingle_well, well_inc)
    Vright_well = np.roll(Vsingle_well, -well_inc)
    Vdouble_well = Vleft_well + Vright_well

    Hdouble_well = KE + 2*V0*np.diag( Vdouble_well )
    eps_double_well, psi_double_well = sorted_eig( Hdouble_well )
    
    R = 2*well_inc*dx
    Rvals[i] = R
    Evals[i, :] = eps_double_well[0:2]
    
    well_inc += 1 # push the wells one step appart 

    # plot everything 

    ax1.cla() # Clear the current axes
    ax2.cla() # Clear the current axes

    ax1.grid()

    colors = ['tab:blue', 'tab:orange']
    for level in [0, 1]:

        wf = np.real( psi_double_well[:, level] )
        ax1.plot( xvals, np.sign( wf[0] )*wf, c = colors[level], label = 'level ' + str(level) + ' (exact)' )


    ax1.set_xlim([ min(xvals), max(xvals) ])
    ax1.set_ylim([ -0.2, 0.2 ])
    ax1.set_title( r'$R = ' + str(np.round( R, 3 )) + '~ a_0$' )
    ax1.set_ylabel( r'The wavefunction' )
    ax1.legend()


    # second plot

    ax2.grid()
    ax2.plot( xvals, Vdouble_well, 'black' )

    for j in [0, 1]:
        ax2.axhline( y = eps_double_well[j], c = colors[j], label = 'level ' + str(level) + ' (exact)' )
    
    ax2.set_ylim([ -1.2, 0.2 ])
    ax2.set_ylabel( r'Energy ($Ry = 13.6$ eV)' )
    ax2.set_xlabel( r'x (Borh radii $a_0$)' )
    plt.tight_layout()

    plt.show()
    plt.pause( 0.05 )


plt.ioff()
fig2 = plt.figure(2, figsize = (8, 4.5))
ax = fig2.add_subplot()
ax.grid( zorder = 1 )

ax.plot( Rvals, Evals[:, 0], color = 'tab:blue', label = 'level 0 exact' )
ax.plot( Rvals, Evals[:, 1], color = 'tab:orange', label = 'level 1 exact' )

ax.legend()
ax.set_xlabel( r'$R$ in Bohr radii $a_0$' )
ax.set_ylabel( r'Energy ($Ry = 13.6$ eV)' )
ax.set_xlim([ 0.0, max(Rvals) ])

plt.show()
