import numpy as np
from test_data import (
    DAE_11,
    Test_0_inputs, Test_0_Results,
    Test_1_inputs, Test_1_Results,
    Test_2_inputs, Test_2_Results,
    Test_5_inputs, Test_5_Results,
    Test_6_inputs, Test_6_Results,
    Test_7_inputs, Test_7_Results,
)
import xfoil as xf


# =============================================================================
#  Helper functions
# =============================================================================

def setup_naca(desig, n_panels=160):
    xf.initialise()
    xf.quiet(True)
    xf.naca_load(desig)
    xf.repanel(n_panels, 1.0, 0.15, 0.2, 1.0, 1.0, 1.0, 1.0)


def setup_airfoil(x, y, n_panels=160):
    xf.initialise()
    xf.quiet(True)
    xf.airfoil(x, y)
    xf.repanel(n_panels, 1.0, 0.15, 0.2, 1.0, 1.0, 1.0, 1.0)


def apply_flow_conditions(inputs):
    xf.setcon(inputs['Re'], inputs['Ma'])
    xf.setiter(300)
    xf.trpars(inputs['XTR Top'], inputs['XTR Bot'], inputs['Ncrit'])


def get_f2py_geometry(n_buf=500):
    x_arr, y_arr, n_out = xf.get_xy(n_buf)
    cp_arr, n_out_cp    = xf.get_cp(n_buf)
    return x_arr[:n_out], y_arr[:n_out], cp_arr[:n_out_cp]


def get_f2py_bl(n_buf=600):
    """Return boundary layer arrays sliced to the actual point count."""
    s, x, y, ue, dstr, thet, cf, hk, cdis, ct, n_out = xf.get_bl(n_buf)
    sl = slice(n_out)
    return (s[sl], x[sl], y[sl], ue[sl], dstr[sl],
            thet[sl], cf[sl], hk[sl], cdis[sl], ct[sl], n_out)


# =============================================================================
#  Test runners — silent, return result dicts
# =============================================================================

def run_single_alpha_test(title, inputs, results, airfoil_type,
                          airfoil_data=None):
    if airfoil_type == 'naca':
        desig = inputs['Airfoil'].replace('NACA', '').replace('naca', '')
        setup_naca(desig)
    else:
        setup_airfoil(airfoil_data['x'], airfoil_data['y'])

    apply_flow_conditions(inputs)
    cl, cd, cm, xtt, xtb, conv = xf.alpha(inputs['alpha'])
    x_f2py, y_f2py, cp_f2py    = get_f2py_geometry()

    x_ref  = results['x']
    y_ref  = results['y']
    cp_ref = results['cp']

    panel_match = len(x_ref) == len(x_f2py)
    if panel_match:
        x_d  = x_f2py  - x_ref
        y_d  = y_f2py  - y_ref
        cp_d = cp_f2py - cp_ref
        geom = {
            'n_ref': len(x_ref), 'n_f2py': len(x_f2py),
            'panel_match': True,
            'x_max':   float(np.max(np.abs(x_d))),
            'x_mean':  float(np.mean(np.abs(x_d))),
            'x_rms':   float(np.sqrt(np.mean(x_d**2))),
            'y_max':   float(np.max(np.abs(y_d))),
            'y_mean':  float(np.mean(np.abs(y_d))),
            'y_rms':   float(np.sqrt(np.mean(y_d**2))),
            'cp_max':  float(np.max(np.abs(cp_d))),
            'cp_mean': float(np.mean(np.abs(cp_d))),
            'cp_rms':  float(np.sqrt(np.mean(cp_d**2))),
        }
    else:
        geom = {
            'n_ref': len(x_ref), 'n_f2py': len(x_f2py),
            'panel_match': False,
        }

    # --- Boundary layer comparison ---
    def _bl_stats(d):
        a = np.asarray(d)
        return {
            'max':  float(np.max(np.abs(a))),
            'mean': float(np.mean(np.abs(a))),
            'rms':  float(np.sqrt(np.mean(a**2))),
        }

    s_f2py, _, _, ue_f2py, dstr_f2py, thet_f2py, cf_f2py, hk_f2py, _, _, n_bl = \
        get_f2py_bl()
    n_ref_bl = len(results['s'])

    if panel_match and n_bl >= n_ref_bl:
        bl = {
            'available': True,
            'n_out': n_bl, 'n_ref': n_ref_bl,
            's':     _bl_stats(s_f2py[:n_ref_bl]    - results['s']),
            'Ue':    _bl_stats(ue_f2py[:n_ref_bl]   - results['Ue_vinf']),
            'Dstar': _bl_stats(dstr_f2py[:n_ref_bl]  - results['Dstar']),
            'Theta': _bl_stats(thet_f2py[:n_ref_bl]  - results['Theta']),
            'Cf':    _bl_stats(cf_f2py[:n_ref_bl]    - results['Cf']),
            'H':     _bl_stats(hk_f2py[:n_ref_bl]    - results['H']),
        }
    else:
        bl = {'available': False, 'n_out': n_bl, 'n_ref': n_ref_bl}

    return {
        'title': title, 'type': 'single',
        'alpha': inputs['alpha'], 'converged': bool(conv),
        'scalars': {
            'CL':  (results['CL'],      float(cl)),
            'CD':  (results['CD'],      float(cd)),
            'CM':  (results['CM'],      float(cm)),
            'XTT': (results['Xtr Top'], float(xtt)),
            'XTB': (results['Xtr Bot'], float(xtb)),
        },
        'geom': geom,
        'bl':   bl,
    }


def run_sweep_test(title, inputs, results, airfoil_type, airfoil_data=None):
    if airfoil_type == 'naca':
        desig = inputs['Airfoil'].replace('NACA', '').replace('naca', '')
        setup_naca(desig)
    else:
        setup_airfoil(airfoil_data['x'], airfoil_data['y'])

    apply_flow_conditions(inputs)

    alpha_i    = inputs['alpha_i']
    alpha_f    = inputs['alpha_f']
    alpha_step = inputs['alpha_step']
    n_alpha    = int(round((alpha_f - alpha_i) / alpha_step)) + 1

    a_arr, cl_arr, cd_arr, cm_arr, xtu_arr, xtb_arr, conv_arr = \
        xf.aseq(alpha_i, alpha_f, n_alpha)

    cl_d, cd_d, cm_d, xtt_d, xtb_d = [], [], [], [], []
    for i, ref_a in enumerate(results['alphas']):
        idx = int(np.argmin(np.abs(a_arr - ref_a)))
        if abs(a_arr[idx] - ref_a) > 0.01:
            continue
        cl_d.append(float(cl_arr[idx]  - results['CLs'][i]))
        cd_d.append(float(cd_arr[idx]  - results['CDs'][i]))
        cm_d.append(float(cm_arr[idx]  - results['CMs'][i]))
        xtt_d.append(float(xtu_arr[idx] - results['Xtr Tops'][i]))
        xtb_d.append(float(xtb_arr[idx] - results['Xtr Bots'][i]))

    def stats(arr):
        a = np.array(arr)
        return (float(np.max(np.abs(a))),
                float(np.mean(np.abs(a))),
                float(np.sqrt(np.mean(a**2))))

    return {
        'title': title, 'type': 'sweep',
        'n_ref': len(results['alphas']), 'n_matched': len(cl_d),
        'alpha_range': (alpha_i, alpha_f, alpha_step),
        'summary': {
            'CL':      stats(cl_d),
            'CD':      stats(cd_d),
            'CM':      stats(cm_d),
            'Xtr Top': stats(xtt_d),
            'Xtr Bot': stats(xtb_d),
        },
    }


# =============================================================================
#  Run all tests
# =============================================================================

print("Running tests...", flush=True)

r0 = run_single_alpha_test(
    "Test 0 — NACA 0012 | Inviscid | Ma=0.0 | alpha=5 deg",
    Test_0_inputs, Test_0_Results, airfoil_type='naca'
)
print("  Test 0 done", flush=True)

r1 = run_single_alpha_test(
    "Test 1 — NACA 0012 | Re=1,000,000 | Ma=0.0 | alpha=5 deg",
    Test_1_inputs, Test_1_Results, airfoil_type='naca'
)
print("  Test 1 done", flush=True)

r2 = run_single_alpha_test(
    "Test 2 — NACA 0012 | Re=1,000,000 | Ma=0.6 | alpha=5 deg",
    Test_2_inputs, Test_2_Results, airfoil_type='naca'
)
print("  Test 2 done", flush=True)

r5 = run_sweep_test(
    "Test 5 — NACA 0012 | Re=650,000 | Ma=0.0 | alpha sweep -15 to 20 deg",
    Test_5_inputs, Test_5_Results, airfoil_type='naca'
)
print("  Test 5 done", flush=True)

r6 = run_single_alpha_test(
    "Test 6 — DAE-11 | Re=250,000 | Ma=0.0 | alpha=5 deg",
    Test_6_inputs, Test_6_Results, airfoil_type='custom', airfoil_data=DAE_11
)
print("  Test 6 done", flush=True)

r7 = run_sweep_test(
    "Test 7 — DAE-11 | Re=450,000 | Ma=0.0 | alpha sweep -10 to 15 deg",
    Test_7_inputs, Test_7_Results, airfoil_type='custom', airfoil_data=DAE_11
)
print("  Test 7 done", flush=True)

print("\nAll tests complete. Printing results...\n")


# =============================================================================
#  Print all results
# =============================================================================

def print_single_result(r):
    print(f"\n{'='*65}")
    print(f"  {r['title']}")
    print(f"{'='*65}")
    print(f"  alpha = {r['alpha']} deg   Converged: {r['converged']}")
    print(f"\n  {'Parameter':<12} {'XFOIL':>12} {'f2py':>12} {'Difference':>12}")
    print("  " + "-" * 52)
    for param, (ref, got) in r['scalars'].items():
        print(f"  {param:<12} {ref:>12.6f} {got:>12.6f} {got - ref:>+12.6f}")

    g = r['geom']
    print(f"\n  Panel count — reference: {g['n_ref']}, f2py: {g['n_f2py']}")
    if not g['panel_match']:
        print("  *** Panel counts differ — paneling mismatch ***")
        return
    print(f"\n  {'Metric':<20} {'x':>12} {'y':>12} {'Cp':>12}")
    print("  " + "-" * 60)
    print(f"  {'Max abs diff':<20} {g['x_max']:>12.6f}"
          f" {g['y_max']:>12.6f} {g['cp_max']:>12.6f}")
    print(f"  {'Mean abs diff':<20} {g['x_mean']:>12.6f}"
          f" {g['y_mean']:>12.6f} {g['cp_mean']:>12.6f}")
    print(f"  {'RMS diff':<20} {g['x_rms']:>12.6f}"
          f" {g['y_rms']:>12.6f} {g['cp_rms']:>12.6f}")

    bl = r['bl']
    print(f"\n  Boundary layer  (n_out={bl['n_out']}, n_ref={bl['n_ref']})")
    if not bl['available']:
        print("  *** BL comparison unavailable — count mismatch ***")
        return
    print(f"  {'Quantity':<10} {'Max|diff|':>12} {'Mean|diff|':>12} {'RMS':>12}")
    print("  " + "-" * 50)
    for name in ('s', 'Ue', 'Dstar', 'Theta', 'Cf', 'H'):
        st = bl[name]
        print(f"  {name:<10} {st['max']:>12.6f} {st['mean']:>12.6f}"
              f" {st['rms']:>12.6f}")


def print_sweep_result(r):
    print(f"\n{'='*65}")
    print(f"  {r['title']}")
    print(f"{'='*65}")
    ai, af, step = r['alpha_range']
    print(f"  Sweep: {ai} to {af} deg, step {step}  "
          f"— matched {r['n_matched']} / {r['n_ref']} reference points")
    print(f"\n  {'Quantity':<12} {'Max|diff|':>12} {'Mean|diff|':>12} {'RMS':>12}")
    print("  " + "-" * 52)
    for name, (mx, mn, rms) in r['summary'].items():
        print(f"  {name:<12} {mx:>12.6f} {mn:>12.6f} {rms:>12.6f}")


print_single_result(r0)
print_single_result(r1)
print_single_result(r2)
print_sweep_result(r5)
print_single_result(r6)
print_sweep_result(r7)