# MATLAB: helium cylinder with parallel valves

## Class

`addition_gas_cylinder_parallel_valve.m` — orifice + two parallel solenoid valves, cylinder blowdown state update via REFPROP.

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
