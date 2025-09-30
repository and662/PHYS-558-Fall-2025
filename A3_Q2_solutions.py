# 558 Question 2, Andrew Mckenna 

import numpy as np 
from matplotlib import pyplot as plt
# from matplotlib import rc
# rc( 'text', usetex = True ) # use nice LaTeX font
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

Hsingle_well = KE + 2*V0*np.diag(Vsingle_well)
eps_single_well, psi_single_well = sorted_eig( Hsingle_well )


# initialize the figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize = (8, 6))


for i in range(Rres):

    Vdouble_well = np.zeros((xres, ))
    
    Vleft_well = np.roll(Vsingle_well, well_inc)
    Vright_well = np.roll(Vsingle_well, -well_inc)
    Vdouble_well = Vleft_well + Vright_well

    Hdouble_well = KE + 2*V0*np.diag( Vdouble_well )
    Hleft_well = KE + 2*V0*np.diag( Vleft_well )
    Hright_well = KE + 2*V0**np.diag( Vright_well )

    eps_double_well, psi_double_well = sorted_eig( Hdouble_well )
    
    eps_single_well_left, psi_single_well_left = sorted_eig( Hleft_well )
    eps_single_well_right, psi_single_well_right = sorted_eig( Hright_well )

    R = 2*well_inc*dx
    Rvals[i] = R
    Evals[i, :] = eps_double_well[0:2]
    

    # pick the lowest energy levels as a basis:
    
    phi_left = np.roll(psi_single_well[:, 0], -well_inc)
    phi_right = np.roll(psi_single_well[:, 0], well_inc)
    phi = np.array([phi_left, phi_right])

    # orthogonalize the basis with Lowdin:
    
    A = phi @ hc( phi )
    eigenvalues, eigenvectors = np.linalg.eig(np.transpose( A )) # recall A is unitary so the SVD = the diagonalization 
    B = eigenvectors @ np.linalg.inv(np.sqrt(np.diag( eigenvalues ))) @ hc(eigenvectors)

    varphi1 = B[0, 0]*phi_left + B[0, 1]*phi_right
    varphi2 = B[1, 0]*phi_left + B[1, 1]*phi_right
    
    # renormalize:
    
    varphi1 = (1.0/np.linalg.norm(varphi1))*varphi1
    varphi2 = (1.0/np.linalg.norm(varphi2))*varphi2
    varphi = np.array([varphi1, varphi2])
    
    # compute the effective Hamiltonian: 

    effective_Ham = varphi @ Hdouble_well @ hc( varphi )
    TBA_energies, TBA_coeff = sorted_eig( effective_Ham )
    Evals_TBA[i, :] = TBA_energies[:]

    TBA_states = np.zeros((xres, 2))
    TBA_states0 = TBA_coeff[0, 0]*varphi1 + TBA_coeff[0, 1]*varphi2
    TBA_states1 = TBA_coeff[1, 0]*varphi1 + TBA_coeff[1, 1]*varphi2


    well_inc += 1 # push the wells one step appart 


    # plot everything 

    ax1.cla() # Clear the current axes
    ax2.cla() # Clear the current axes

    ax1.grid()

    colors = ['tab:blue', 'tab:orange']
    for level in [0, 1]:

        wf = np.real( psi_double_well[:, level] )
        ax1.plot( xvals, np.sign( wf[0] )*wf, c = colors[level], label = f'level { level } (exact)'  )


    # plot the basis

    ax1.plot( xvals, np.sign( TBA_states1[0] )*TBA_states1, 'tab:green', label = 'level 0 (TBA)' )
    ax1.plot( xvals, np.sign( TBA_states0[0] )*TBA_states0, 'tab:red', label = 'level 1 (TBA)' )
    
    ax1.set_xlim([ min(xvals), max(xvals) ])
    ax1.set_ylim([ -0.2, 0.2 ])
    ax1.set_title( r'$R = ' + f'{ np.round( R, 3 ) }' '~ a_0$' )
    ax1.set_ylabel( r'The wavefunction' )
    ax1.legend()


    # second plot

    ax2.grid()
    ax2.plot( xvals, Vdouble_well, 'black' )

    ax2.axhline( y = eps_single_well[0], color = 'grey', linestyle = '--', label = 'single well' )
    ax2.axhline( y = TBA_energies[0], color = 'tab:red', label = 'single well' )
    ax2.axhline( y = TBA_energies[1], color = 'tab:green', label = 'single well' )
    
    for j in [0, 1]:
        ax2.axhline( y = eps_double_well[j], c = colors[j], label = f'level { level } (exact)' )
    
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

ax.plot( Rvals, Evals_TBA[:, 0], color = 'tab:red', label = 'level 0 TBA' )
ax.plot( Rvals, Evals_TBA[:, 1], color = 'tab:green', label = 'level 1 TBA' )

ax.axhline( y = eps_single_well[0], color = 'grey', linestyle = '--', label = 'single well' )

ax.legend()
ax.set_xlabel( r'$R$ in Bohr radii $a_0$' )
ax.set_ylabel( r'Energy ($Ry = 13.6$ eV)' )
ax.set_xlim([ 0.0, max(Rvals) ])

plt.show()
