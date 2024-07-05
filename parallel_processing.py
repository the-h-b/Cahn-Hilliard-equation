from petsc4py import PETSc
from mpi4py import MPI
from dolfin import *
import numpy as np

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Define mesh and function space
nx, ny = 50, 50
mesh = UnitSquareMesh(comm, nx, ny)
degree = 1
P1 = FiniteElement('Lagrange', triangle, degree)
ME = FunctionSpace(mesh, MixedElement([P1, P1]))

# Define trial and test functions
u = Function(ME)
v = TestFunction(ME)
c, mu = split(u)
vc, vmu = split(v)

# Define parameters
epsilon = 5e-2
M = 1.0

# Define initial conditions
class InitialConditions(UserExpression):
    def eval(self, values, x):
        values[0] = 0.63 + 0.02 * (0.5 - np.random.rand())
        values[1] = 0.0

    def value_shape(self):
        return (2,)

u.interpolate(InitialConditions(degree=degree))

# Define the Cahn-Hilliard free energy and chemical potential
c = variable(c)
f = (1 / 4) * (c*2 - 1)*2
dfdc = diff(f, c)

# Weak statement of the equations
L0 = c * vmu * dx - M * dot(grad(mu), grad(vmu)) * dx
L1 = mu * vc * dx - (dfdc - epsilon**2 * div(grad(c))) * vc * dx
L = L0 + L1

# Derivative of F
a = derivative(L, u)

# Setup nonlinear solver
problem = NonlinearVariationalProblem(L, u, J=a)
solver = NonlinearVariationalSolver(problem)

# Solver parameters
prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1e-6
prm['newton_solver']['relative_tolerance'] = 1e-5
prm['newton_solver']['maximum_iterations'] = 100
prm['newton_solver']['relaxation_parameter'] = 1.0

# Specify a different linear solver and preconditioner
prm['newton_solver']['linear_solver'] = 'gmres'
prm['newton_solver']['preconditioner'] = 'ilu'

# Enable more verbose solver output
prm['newton_solver']['report'] = True
prm['newton_solver']['error_on_nonconvergence'] = False

# Time-stepping
T = 5.0  # final time
num_steps = 50  # number of time steps
dt = T / num_steps  # time step size
t = 0

# Time-stepping loop
for n in range(num_steps):
    # Update current time
    t += dt
    if rank == 0:
        print(f'Time step {n}, time {t}')

    # Solve the problem with error handling
    try:
        solver.solve()
    except RuntimeError as e:
        if rank == 0:
            print(f"Error in solving at time step {n}: {e}")
        break

    # Extract solutions
    (c, mu) = u.split()

    # Print information for debugging
    if rank == 0:
        print(f"Norm of c: {c.vector().norm('l2')}, Norm of mu: {mu.vector().norm('l2')}")

# Save solution in VTK format
if rank == 0:
    file = File('cahn_hilliard_solution.pvd')
    file << u