# MATLAB: helium cylinder with parallel valves

## Classes

- `addition_gas_cylinder_parallel_valve.m` — orifice + two parallel solenoid valves, cylinder blowdown via REFPROP + wall heat transfer.
- `gas_cylinder_model.m` — thin 3-arg alias for scripts that call `gas_cylinder_model(ic, orifice, [v1,v2])`.

### Construct (matches user script)

```matlab
gas_initial_conditions = struct();
gas_initial_conditions.P0_init = 35e6;  % Pa (alias of PO_init)
gas_initial_conditions.T0_init = 300;   % K (alias of TO_init)
gas_initial_conditions.V = 0.18;        % m^3

orifice = struct('Cd', 1, 'A', pi*(0.005)^2);   % as written: 5 mm radius
valve1 = struct('Cd', 0.9, 'A', pi*(0.001)^2);
valve2 = struct('Cd', 0.9, 'A', pi*(0.001)^2);

cyl = gas_cylinder_model(gas_initial_conditions, orifice, [valve1, valve2]);
```

`P0_init`/`T0_init` and `PO_init`/`TO_init` are both accepted in `mass_flow`. See `demo_gas_cylinder_construct.m`.

Add the `matlab/` folder to the MATLAB path (ahead of any older `gas_cylinder_model`) so the 3-arg alias is found.

## Physics check: `test_cylinder_blowdown.m`

Runnable blowdown test with the user init (35 MPa, 300 K, 0.18 m³, orifice + two valves):

```matlab
addpath('<repo>/matlab');   % ahead of any older gas_cylinder_model
test_cylinder_blowdown
```

**What it does**

1. Constructs `gas_cylinder_model` and smokes `networkFlow` (no REFPROP required).
2. If `refpropm` is available: steps `mass_flow` until \(P \approx 0.5\,\mathrm{MPa}\), logs \(P,T,\dot m\) vs time into workspace `blowdown_log`.
3. Asserts **pass** when final \(T > 200\,\mathrm{K}\) (not tens of kelvin) and prints comparison to an isentropic estimate (ideal-gas and, if possible, REFPROP \(S=\mathrm{const}\)).
4. If REFPROP is missing: prints a clear **SKIP** for the state-update check; construct + `networkFlow` still run.

**Expected with defaults + REFPROP:** \(T\) at 0.5 MPa stays well above 200 K (often nearer ambient), while isentropic estimates are ~tens of K.

## Physics change (isentropic → energy + wall heat transfer)

The previous cylinder update forced **constant entropy** after mass left:

```matlab
TO_new = fzero(@(T) refpropm('S',...) - entropy0, TO);
```

Pure isentropic expansion of high-pressure helium down to **0.5 MPa** drives gas temperature to **tens of kelvin**, which is unphysical for a real bottle that exchanges heat with the wall and ambient.

### New control-volume update

Over each `time_step`:

\[
U_{n+1} = m_n u_n - \Delta m\, h_{\mathrm{exit}} + \dot Q\,\Delta t
\]

with \(\dot Q = h_{\mathrm{conv}} A_{\mathrm{wall}} (T_{\mathrm{wall}} - T_{\mathrm{gas}})\), exit enthalpy from the bottle state, and \(T_{n+1}\) from REFPROP \(U(T,\rho_{n+1})\).

Optional lumped wall:

\[
m_{\mathrm{wall}} c_{\mathrm{wall}} \frac{d T_{\mathrm{wall}}}{dt}
  = -\dot Q + h_{\mathrm{amb}} A_{\mathrm{ext}} (T_{\mathrm{amb}} - T_{\mathrm{wall}})
\]

Default `m_wall = 0` keeps an **isothermal wall at `T_amb`** (strong thermal coupling limit).

Nozzle `T_out` in `elementFlow` remains isentropic (orifice/valve path only); it does not set bottle temperature.

## Default heat-transfer parameters

| Parameter | Default | Units | Meaning |
|-----------|---------|-------|---------|
| `h_conv` | `100` | W/(m²·K) | Gas–inner-wall convection |
| `A_wall` | from `V` (L/D≈4) | m² | Inner wall area |
| `T_amb` | `293.15` | K | Ambient |
| `T_wall` | `T_amb` | K | Wall temperature |
| `m_wall` | `0` | kg | `0` ⇒ isothermal wall at `T_amb` |
| `c_wall` | `500` | J/(kg·K) | Steel-like wall specific heat |
| `h_amb` | `10` | W/(m²·K) | Outer wall–ambient convection |
| `A_ext` | `A_wall` | m² | Outer area |

Override via fields on `initial_conditions` (and `T_wall` is carried in `new_conditions`).

## Expected temperature at 0.5 MPa

- **Old isentropic:** tens of K from a typical high-P / near-ambient start.
- **New defaults (isothermal wall ≈ ambient, moderate `h_conv`):** for slow-to-moderate discharge, bottle gas stays **near ambient** (often ~250–293 K rather than cryogenic). Faster discharge or weaker `h_conv` cools more; raise `h_conv` / `A_wall` or keep `m_wall = 0` to stay closer to isothermal.
