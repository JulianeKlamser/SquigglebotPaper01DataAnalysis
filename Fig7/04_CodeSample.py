import numpy as np
import matplotlib.pyplot as plt
import argparse
import time

def read_command_line_args():
    parser = argparse.ArgumentParser(description="Read command-line arguments separated by spaces.")
    parser.add_argument("args", nargs="+", help="List of arguments separated by spaces.")
    return parser.parse_args()
    
def wca_forceAmplitude(r, sigma_ij, cutOffScale):
    """Compute WCA force for interparticle distances r (can be a NumPy array)."""
    r = np.asarray(r)  # Ensure r is a NumPy array
    mask = (r < (cutOffScale * sigma_ij)) & (r > 0)
    force = np.zeros_like(r)
    force[mask] = 4.0 * (12 * (sigma_ij[mask]**12) / r[mask]**14 - 6 * (sigma_ij[mask]**6) / r[mask]**8)
    return force
        
def SafeConf(Path, N, dt, t_tot, Conf_ID, num_steps, positions, velocities, theta, radii, CPUTime):
    DataFormat = np.zeros((N,6))
    DataFormat[:,:2] = positions
    DataFormat[:,2:4] = velocities
    DataFormat[:,4] = theta
    DataFormat[:,5] = radii
    
    ConfFile = Path + '/Particles%d.txt' % Conf_ID
    np.savetxt( ConfFile, DataFormat, delimiter=' ', newline='\n', header=' x y vx vy theta radius')
    
    ConfFile = Path + '/Particles.txt'
    np.savetxt( ConfFile, DataFormat, delimiter=' ', newline='\n', header=' x y vx vy theta radius')
    
    TimeFile = Path + '/Time.txt'
    with open(TimeFile, "a") as myfile:
        myfile.write("%d %.8f %d %.8f %.8f %s\n" % (Conf_ID, t_tot, num_steps, dt, CPUTime, "hostName"))
    

# Parameters (dummies - will be read later)
N = 30                  # Number of particles
L = 10.0                # Box size (periodic boundaries)
t_tot = 0.0             # total simulation time at start of simulation
Conf_ID = 0             # total number of already saved configurations
dt = 0.005              # Time step
Delta_t = 10            # simulation time between configuration snapshots
N_conf = 100            # number of configuration snapshots
D_r = 0.1               # Rotational diffusion
v_0 = 1.0               # Self-propulsion speed
#epsilon = 1.0           # WCA potential strength
cutOffScale = 2.**(1./6.)
DisplacementCutOff = 0.15
Path = '/Users/jklamser/Documents/Ludovic/LinearAlignmentCppCode/Data/N000030L7.091678/v_o1.000000D_r0.100000/gamma0.000000alignCutOffBySigma2.000000'

# Initialize particle positions and orientations (also dummies)
positions = np.random.uniform(0, L, (N, 2))
theta = np.random.uniform(0, 2*np.pi, N)
sigma = np.ones(N)

# read command line arguments
print('\n\n--------\nrun code as\npython3 Code.py 100 /path/to/dat/folder/with/init/files\nwhere 100 is the number of configuraitons to safe\n--------\n')
ARGS = read_command_line_args()
assert(len(ARGS.args) == 2)
N_conf = int(ARGS.args[0])
Path = ARGS.args[1]
''' Note that there should be a verification of the command line arguments here '''

# Read time file
TimeFile = Path + '/Time.txt'
Data = np.loadtxt(TimeFile, usecols=(0,1), skiprows=1)
Conf_ID, t_tot = Data[-1,0], Data[-1,1]

# read Parameters
ParamFile = Path + '/Parameters.txt'
N, dummy1, dummy2, dt, Delta_t, f_a, D_r, gamma = np.loadtxt(ParamFile)
N = int(N)

# read periods
PeriodsFile = Path + '/Periods.txt'
Periods = np.loadtxt(PeriodsFile)
np.testing.assert_array_equal(Periods[0], Periods[1])
assert (Periods[0,0] == 0)
L = Periods[0,1]
RadiusOfWall = (L - 5.) / 2.
CircularWallCenter = np.array([RadiusOfWall, RadiusOfWall])

# read pos
ConfFile = Path + '/Particles.txt'
Data = np.loadtxt(ConfFile)
positions = Data[:,:2]
velocities = Data[:,2:4]
theta = Data[:,4]
radii = Data[:,5]
sigma = 2*radii
DisplacementCutOff = 0.15*min(sigma)

# Simulation loop
num_steps = int(Delta_t / dt)
for conf_i in range(N_conf):
    CPUTime_Start = time.process_time()
    for step_j in range(num_steps):
        # init forces to zero
        forces = np.zeros((N, 2))
        
        # Compute interaction forces
        for i in range(N):
            rij = positions[i] - positions
            # rij -= L * np.round(rij / L)  # NO periodic boundary conditions needed
            r = np.linalg.norm(rij, axis = 1)
            sigma_ij = radii[i] + radii
            F = wca_forceAmplitude(r, sigma_ij, cutOffScale)[:, np.newaxis] * rij # force on i
            forces[i] += np.sum(F, axis = 0)
        
        # get force of wall on particles
        r_center_vec = CircularWallCenter - positions # points towards center
        centerDist = np.linalg.norm(r_center_vec, axis = 1)
        wallDist = RadiusOfWall - centerDist
        forces += wca_forceAmplitude(wallDist, radii, cutOffScale)[:, np.newaxis] * r_center_vec # force on i
        
        # add the active force on particle 0
        forces[0] += f_a * np.array([np.cos(theta[0]), np.sin(theta[0])]) # only particle 0 has non-zero self-propulsion velocity
        
        # Update positions, velocities and orientations
        # update positions r(t+dt) = r(t) + dt v(t)
        disp = dt * velocities
        if( ( np.abs(disp) > DisplacementCutOff ).any() ):
            print('displacement to large ! Reduce dt !')
            print('at least one displacement component is larger than %f' % DisplacementCutOff)
            exit()
        positions += disp
        # positions %= L  # NO periodic boundary conditions needed
        
        # Update velocities
        # v(t+dt) = (1 - dt gamma) * v(t) + dt F(t), for active particle F = F_steric + F_active
        velocities = (1. - dt * gamma) * velocities  + dt * forces
        
        # update orientation of active particle
        dtheta = np.sqrt(2 * D_r * dt) * np.random.randn()
        theta[0] = np.mod(theta[0] + dtheta, 2*np.pi)
        
        
    CPUTime_END = time.process_time()
    CPUTime = (CPUTime_END - CPUTime_Start)/60.
    t_tot += dt*num_steps
    Conf_ID += 1
    # safe configuration
    SafeConf(Path, N, dt, t_tot, Conf_ID, num_steps, positions, velocities, theta, radii, CPUTime)
        

