


# Hur får man ut tidsderivatan av ordning 1 och 2?
# i.e u_t och u_tt



# Diffusion term
dterm = d * dt/2 *(Tdx * unew + Tdx * uold)

# The LW-scheme
# u(t + dt, x) = u(t,x) + dt*u_t + dt^2 / 2! * u_tt 
def LW(uold):
    unew = uold + dt * u_t + dt**2 / 2 * u_tt
    return unew


# Simple second order discretization of the viscous Burgers equation
# The nonlinear convection part is treated by the explicit LW method and the linear diffusion part is taken care of by the 
# trapezoidal Rule

def scheme4(uold):
    A = LaxWen(uold) + d * dt/2 * Tdx * uold
    B = (np.identity(len(uold)) - d * dt/2 * Tdx)
    B_inv = np.linalg.inv(B)
    unew = np.dot(B_inv, A)
    