from dolfin import *

# Create mesh and define function space
mesh = RectangleMesh(Point(0.0, 0.0), Point(2.0, 1.0), 64, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundaries
def left_boundary(x):
    return x[0] < DOLFIN_EPS


def right_boundary(x):
    return x[0] > 2.0 - DOLFIN_EPS


# Define boundary conditions
u_left = Constant(0.0)
u_right = Constant(1.0)

bc_left = DirichletBC(V, u_left, left_boundary)
bc_right = DirichletBC(V, u_right, right_boundary)
bcs = [bc_left, bc_right]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

# Two heat sources inside the plate
f = Expression(
    "25*exp(-(pow(x[0] - 0.6, 2) + pow(x[1] - 0.3, 2)) / 0.02)"
    " + 30*exp(-(pow(x[0] - 1.4, 2) + pow(x[1] - 0.7, 2)) / 0.03)",
    degree=2
)

a = inner(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Plot solution and save image
import matplotlib.pyplot as plt

image = plot(u)
plt.colorbar(image, label="Temperature")
plt.title("Temperature in a plate with two heat sources")
plt.savefig("poisson_temperature.png", dpi=200, bbox_inches="tight")
plt.show()