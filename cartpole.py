import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import solve_discrete_are
from scipy.optimize import minimize
from scipy.signal import cont2discrete

# ---- 1. Model: cart mass M, pendulum point mass m on a rod of length l ----
M, m, l, g = 1.0, 0.2, 0.5, 9.81
dt = 0.02                       # 50 Hz control loop
U_MAX = 10.0                    # actuator limit (N)

# Linearized around upright, state x = [cart pos, cart vel, angle, angular vel]
A = np.array([[0, 1, 0, 0],
              [0, 0, -m*g/M, 0],
              [0, 0, 0, 1],
              [0, 0, (M+m)*g/(M*l), 0]])
B = np.array([[0], [1/M], [0], [-1/(M*l)]])
Ad, Bd, *_ = cont2discrete((A, B, np.eye(4), np.zeros((4, 1))), dt)

def plant(x, u):
    """Nonlinear cart-pole, one RK4 step (the 'real' system)."""
    def f(s):
        _, xd, th, thd = s
        sin, cos = np.sin(th), np.cos(th)
        den = M + m*sin**2
        xdd = (u + m*sin*(l*thd**2 - g*cos)) / den
        thdd = (-u*cos - m*l*thd**2*sin*cos + (M+m)*g*sin) / (l*den)
        return np.array([xd, xdd, thd, thdd])
    k1 = f(x); k2 = f(x + dt/2*k1); k3 = f(x + dt/2*k2); k4 = f(x + dt*k3)
    return x + dt/6*(k1 + 2*k2 + 2*k3 + k4)

# ---- 2. LQR ----
Q = np.diag([10, 1, 100, 1]); R = np.array([[0.1]])
P = solve_discrete_are(Ad, Bd, Q, R)
K = np.linalg.solve(R + Bd.T @ P @ Bd, Bd.T @ P @ Ad)
def lqr(xhat, u_max=U_MAX):
    return float(np.clip((-K @ xhat)[0], -u_max, u_max))

# ---- 3. MPC with an explicit force limit (condensed QP, solved with SciPy) ----
N = 30                                   # horizon: 30 steps = 0.6 s
Sx = np.vstack([np.linalg.matrix_power(Ad, k) for k in range(1, N+1)])      # (4N x 4)
Su = np.zeros((4*N, N))                                                    # (4N x N)
for k in range(N):
    for j in range(k+1):
        Su[4*k:4*k+4, j] = (np.linalg.matrix_power(Ad, k-j) @ Bd).ravel()
Qbar = np.kron(np.eye(N), Q); Qbar[-4:, -4:] = P   # terminal cost = LQR's P
H = Su.T @ Qbar @ Su + np.kron(np.eye(N), R)
u_prev = np.zeros(N)
def mpc(xhat, u_max=U_MAX):
    global u_prev
    f = Su.T @ Qbar @ (Sx @ xhat)
    res = minimize(lambda U: U @ H @ U + 2*f @ U, u_prev, jac=lambda U: 2*(H @ U + f),
                   bounds=[(-u_max, u_max)]*N, method="L-BFGS-B")
    u_prev = np.r_[res.x[1:], res.x[-1]]          # warm start next step
    return float(res.x[0])

# ---- 4. Kalman filter: we only measure cart position and angle, with noise ----
C = np.array([[1, 0, 0, 0], [0, 0, 1, 0]])
Rn = np.diag([0.01**2, np.deg2rad(1.0)**2])   # sensor noise: 1 cm, 1 degree
Qn = 1e-4 * np.eye(4)                          # process noise
def kf_step(xhat, Pk, u, y):
    xhat = Ad @ xhat + (Bd * u).ravel(); Pk = Ad @ Pk @ Ad.T + Qn
    S = C @ Pk @ C.T + Rn; Kk = Pk @ C.T @ np.linalg.inv(S)
    xhat = xhat + Kk @ (y - C @ xhat); Pk = (np.eye(4) - Kk @ C) @ Pk
    return xhat, Pk

# ---- 5. Run both controllers from the same start (pendulum tilted 20 degrees) ----
def run(ctrl, u_max=U_MAX, start_deg=20, T=5.0, seed=0):
    global u_prev
    u_prev = np.zeros(N)
    rng = np.random.default_rng(seed)
    x = np.array([0, 0, np.deg2rad(start_deg), 0])
    y0 = C @ x + rng.multivariate_normal(np.zeros(2), Rn)
    xhat = np.array([y0[0], 0, y0[1], 0]); Pk = 0.1*np.eye(4)   # start the filter from the first reading
    log = []
    for i in range(int(T/dt)):
        y = C @ x + rng.multivariate_normal(np.zeros(2), Rn)   # noisy sensors
        u = ctrl(xhat, u_max)
        xhat, Pk = kf_step(xhat, Pk, u, y)
        log.append([i*dt, *x, u, *xhat]); x = plant(x, u)
    return np.array(log)    # columns: t, x, xdot, theta, thetadot, u, 4 estimates

def settling_time(log, tol_deg=3.0):
    th = np.rad2deg(np.abs(log[:, 3])); outside = np.where(th > tol_deg)[0]
    if len(outside) == 0: return 0.0
    if outside[-1] == len(log) - 1: return float("inf")   # never settled
    return log[outside[-1]+1, 0]

def balanced(log):
    return np.rad2deg(np.abs(log[-50:, 3])).max() < 5    # still upright in the last second

results = {name: run(c) for name, c in [("LQR", lqr), ("MPC", mpc)]}
for name, log in results.items():
    err = np.rad2deg(np.sqrt(np.mean((log[:, 8] - log[:, 3])**2)))
    print(f"{name}: settles within 3 deg in {settling_time(log):.2f} s, "
          f"peak force {np.abs(log[:, 5]).max():.1f} N, Kalman angle error (RMS) {err:.2f} deg")

# Force-limit sweep: the weakest motor each controller can still recover with
for name, c in [("LQR", lqr), ("MPC", mpc)]:
    ok = [u for u in np.arange(10.0, 4.5, -0.5) if balanced(run(c, u_max=u))]
    print(f"{name}: recovers from 20 deg with a force limit as low as {min(ok):.1f} N" if ok else f"{name}: never recovers")

fig, ax = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
for name, log in results.items():
    ax[0].plot(log[:, 0], np.rad2deg(log[:, 3]), label=name)
    ax[1].plot(log[:, 0], log[:, 5], label=name)
ax[0].set_ylabel("pendulum angle (deg)"); ax[1].set_ylabel("cart force (N)"); ax[1].set_xlabel("time (s)")
ax[0].legend(); plt.tight_layout(); plt.savefig("lqr_vs_mpc.png", dpi=150)
print("saved lqr_vs_mpc.png")
