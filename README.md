# Cart-Pole: LQR vs MPC with a Kalman Filter

Balances a simulated inverted pendulum on a cart with two controllers (LQR and MPC)
while a Kalman filter estimates the full state from noisy position and angle sensors.

Runs entirely on a Mac (tested; about 5 seconds).

## Step-by-step (Mac M4)

1. Open **Terminal** and go into this folder:
   ```bash
   cd ~/robotics-projects/project2-cartpole-mpc
   ```
2. Create and turn on a virtual environment (do step 2's `source` line every time you open a new Terminal):
   ```bash
   python3.11 -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -r requirements.txt
   ```
3. Run it:
   ```bash
   python cartpole.py
   ```
4. You should see 4 printed lines and a new file `lqr_vs_mpc.png`. Open it with:
   ```bash
   open lqr_vs_mpc.png
   ```

## What the script does

| Part | What it is |
|---|---|
| `plant()` | The "real robot": nonlinear cart-pole simulated with RK4 |
| `A, B` -> `Ad, Bd` | Linearized state-space model, discretized at 50 Hz |
| `lqr()` | Fixed-gain LQR controller, clipped at the force limit |
| `mpc()` | 30-step MPC with the force limit inside the optimization |
| `kf_step()` | Kalman filter: measures only cart position (1 cm noise) and angle (1 deg noise) |
| sweep | Finds the weakest motor each controller can still recover with |

## Make it your own (pick at least one, then report what you found)

1. **Track limit:** `pip install cvxpy`, rewrite `mpc()` with a constraint `|cart position| <= 0.3 m`,
   and show MPC respects it while LQR does not.
2. **Horizon:** change `N` to 10, 30, 60 and compare settling time vs. solve time.
3. **Robustness:** change `m` inside `plant()` only (the controller still assumes 0.2 kg) and find where each fails.
4. **Noise:** increase the numbers in `Rn` and plot how the Kalman filter's error grows.

## Put it on GitHub

```bash
git init && git add . && git commit -m "Cart-pole LQR vs MPC with Kalman filter"
gh repo create cartpole-lqr-mpc --public --source=. --push
```
Add `lqr_vs_mpc.png` and your printed results to this README before pushing.

## Resume line (fill in your real numbers)

Modeled a cart-pole in state space and compared LQR and MPC controllers, with a Kalman filter estimating
full state from noisy position and angle sensors ([X] deg RMS error); [your experiment's finding].
