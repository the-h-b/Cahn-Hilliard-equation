from dolfin import *
import numpy as np
import time

# Define mesh and function space
nx, ny = 50, 50
mesh = UnitSquareMesh(nx, ny)
degree = 4
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

# Setup nonlinear problem and solver
problem = NonlinearVariationalProblem(L, u, bcs=None, J=a)
solver = NonlinearVariationalSolver(problem)

# Parameters for solver
solver.parameters['nonlinear_solver'] = 'newton'
solver.parameters['newton_solver']['relative_tolerance'] = 1e-8
solver.parameters['newton_solver']['maximum_iterations'] = 10

# Time-stepping
T = 5.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size
t = 0

# Measure execution time
start_time = time.time()

# Time-stepping loop
for n in range(num_steps):
    # Update current time
    t += dt
    print(f'Time step {n}, time {t}')
    
    # Solve the problem
    solver.solve()
    
    # Extract solutions
    (c, mu) = u.split()
    
    # Update previous solution
    u.vector()[:] = u.vector()

# Calculate total execution time
end_time = time.time()
execution_time = end_time - start_time
print(f"Total execution time: {execution_time:.2f} seconds")

# Save solution in VTK format
file = File('cahn_hilliard_solution.pvd')
file << u