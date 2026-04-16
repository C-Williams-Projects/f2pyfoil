import numpy as np
import sqlite3
import xfoil as xf

DAE_11 = {
    'x': np.array([1,0.986483,0.970008,0.947875,0.918609,0.887967,0.857951,0.82795,0.79778,0.767351,0.736604,0.705653,0.675223,0.645437,0.616266,0.587581,0.559141,0.530681,0.502047,0.473223,0.444274,0.415255,0.386183,0.357092,0.328056,0.299149,0.270425,0.24194,0.213791,0.186118,0.159095,0.132934,0.107915,0.084408,0.062879,0.043966,0.028364,0.016496,0.008293,0.003292,0.000757,0,0.000492,0.002012,0.004415,0.008086,0.014115,0.024573,0.041451,0.064895,0.092608,0.122584,0.153678,0.185384,0.21745,0.249739,0.282172,0.314611,0.346965,0.37924,0.411434,0.443554,0.475605,0.507575,0.539469,0.571312,0.603112,0.634882,0.666635,0.698376,0.730101,0.761789,0.793421,0.824998,0.856562,0.888329,0.920015,0.947994,0.969651,0.986283,1]),
    'y': np.array([0,0.002538,0.005615,0.010045,0.01638,0.023564,0.031128,0.039152,0.047605,0.056417,0.065502,0.074738,0.08372,0.092205,0.099991,0.106907,0.112838,0.117774,0.121778,0.124888,0.127109,0.128434,0.12885,0.128382,0.127043,0.124802,0.121658,0.117585,0.112558,0.106551,0.09954,0.091527,0.082517,0.072565,0.061804,0.05054,0.0393,0.028727,0.019335,0.011468,0.005143,0,-0.003486,-0.006428,-0.008791,-0.0107,-0.012375,-0.013878,-0.014889,-0.015066,-0.014462,-0.013293,-0.011747,-0.009967,-0.008052,-0.006069,-0.004061,-0.002068,-0.000119,0.001752,0.003536,0.0052,0.006734,0.008115,0.009304,0.010285,0.011045,0.01157,0.011866,0.011931,0.011773,0.011393,0.010774,0.0099,0.008759,0.007331,0.005614,0.003868,0.00235,0.001077,0])
}


def Test(file_db, ID, tol=1e-6):

    con = sqlite3.connect(file_db)
    cur = con.cursor()
 
    # inputs
    row = cur.execute("""
        SELECT airfoil, run_type, re, mach, ncrit_top, ncrit_bot,
               xtrip_top, xtrip_bot, input_start, input_end, n_points
        FROM inputs WHERE test_id = ?
    """, (ID,)).fetchone()
    airfoil, run_type, Re, M, ncrit_top, ncrit_bot, \
        xtrip_top, xtrip_bot, input_start, input_end, n_points = row
 
    # polar
    polar_rows = cur.execute("""
        SELECT alpha, cl, cd, cm, xtt, xtb
        FROM polar WHERE test_id = ? ORDER BY rowid
    """, (ID,)).fetchall()
    polar      = np.array(polar_rows, dtype=float)
    ref_alpha  = polar[:, 0]
    ref_cl     = polar[:, 1]
    ref_cd     = polar[:, 2]
    ref_cm     = polar[:, 3]
    ref_xtt    = polar[:, 4]
    ref_xtb    = polar[:, 5]
 
    ref_cp = ref_cf = ref_hstar = None
    ref_s  = ref_x  = ref_y    = None
    ref_ue = ref_dstar = ref_theta = ref_h = None
 
    if run_type in ("alfa", "cl"):
        cp_rows = cur.execute("""
            SELECT cp FROM bl_surface
            WHERE test_id = ? ORDER BY rowid
        """, (ID,)).fetchall()
        ref_cp = np.array(cp_rows, dtype=float).ravel()
 
        if Re != 0.0:
            bls_rows  = cur.execute("""
                SELECT cf, hstar FROM bl_surface
                WHERE test_id = ? ORDER BY rowid
            """, (ID,)).fetchall()
            bls       = np.array(bls_rows, dtype=float)
            ref_cf    = bls[:, 0]
            ref_hstar = bls[:, 1]
 
            bl_rows   = cur.execute("""
                SELECT s, x, y, ue, dstar, theta, h FROM bl
                WHERE test_id = ? ORDER BY rowid
            """, (ID,)).fetchall()
            bl        = np.array(bl_rows, dtype=float)
            ref_s     = bl[:, 0]
            ref_x     = bl[:, 1]
            ref_y     = bl[:, 2]
            ref_ue    = bl[:, 3]
            ref_dstar = bl[:, 4]
            ref_theta = bl[:, 5]
            ref_h     = bl[:, 6]
 
    con.close()
 
    Description = run_type
 
    xf.init()
    xf.quiet(True)
    if airfoil == 'DAE-11':
        xf.setgeom(DAE_11['x'], DAE_11['y'])
    else:
        naca_number = airfoil.upper().replace('NACA', '').strip()
        xf.setnaca(naca_number)
    xf.panel()
 
    xf.flowcons(Re, M)
    xf.settrip(xtrip_top, xtrip_bot)
    xf.setncrit12(ncrit_top, ncrit_bot)
    xf.maxiter(600)
 
    if Description == 'alfa':
        cl, cd, cm, xtt, xtb, conv = xf.alpha(input_start)
        conv_statment = "Converged" if conv else "Not Converged"
        print("")
        print(f"Test {ID}: {Description} - {conv_statment}")
        print("="*115)
        print(f"{'f2pyfoil Scalars:' :<20}" f"CL: {cl:.10f}, CD: {cd:.10f}, CM: {cm:.10f}, Xtc Top: {xtt:.10f}, Xtc Bot: {xtb:.10f}")
        print(f"{'Difference:' :<20}" f"CL: {cl - ref_cl[0]:10.5e}, CD: {cd - ref_cd[0]:10.5e}, CM: {cm - ref_cm[0]:10.5e}, Xtc Top: {xtt - ref_xtt[0]:10.5e}, Xtc Bot: {xtb - ref_xtb[0]:10.5e}")
 
        N = xf.getn(); NBL = xf.getnbl()
        cp = xf.getcp(N); cf = xf.getcf(N)
        s, x, y, ue, dstr, theta, tstar, h, hstari = xf.getmorebl(NBL)
        hstar = hstari[: N]
 
        if Re == 0.0:
            print("Cp Data Summary:")
            print(f"{'':>18}{'cp':>12}")
            print(f"{'Max Difference:':<18} {np.max(np.abs(cp - ref_cp)):>12.5e}")
            print(f"{'Mean Difference:':<18} {np.mean(np.abs(cp - ref_cp)):>12.5e}")
            print(f"{'Tolerance Check:':<18} {str(np.allclose(cp, ref_cp, atol=tol)):>12}")
        else:
            print("BL Data Summary:")
            print("="*115)
            print(f"{'':>18} {'s':>12} {'x':>12} {'y':>12} {'Cp':>12} {'Ue/Vinf':>12} {'dstr':>12} {'theta':>12} {'Cf':>12} {'H':>12} {'H*':>12}")
            print(f"{'Max Difference:':<18} {np.max(np.abs(s - ref_s)):>12.5e} {np.max(np.abs(x - ref_x)):>12.5e} {np.max(np.abs(y - ref_y)):>12.5e} {np.max(np.abs(cp - ref_cp)):>12.5e} {np.max(np.abs(ue - ref_ue)):>12.5e} {np.max(np.abs(dstr - ref_dstar)):>12.5e} {np.max(np.abs(theta - ref_theta)):>12.5e} {np.max(np.abs(cf - ref_cf)):>12.5e} {np.max(np.abs(h - ref_h)):>12.5e} {np.max(np.abs(hstar - ref_hstar)):>12.5e}")
            print(f"{'Mean Difference:':<18} {np.mean(np.abs(s - ref_s)):>12.5e} {np.mean(np.abs(x - ref_x)):>12.5e} {np.mean(np.abs(y - ref_y)):>12.5e} {np.mean(np.abs(cp - ref_cp)):>12.5e} {np.mean(np.abs(ue - ref_ue)):>12.5e} {np.mean(np.abs(dstr - ref_dstar)):>12.5e} {np.mean(np.abs(theta - ref_theta)):>12.5e} {np.mean(np.abs(cf - ref_cf)):>12.5e} {np.mean(np.abs(h - ref_h)):>12.5e} {np.mean(np.abs(hstar - ref_hstar)):>12.5e}")
            print(f"{'Tolerance Check:':<18} {str(np.allclose(s, ref_s, atol=tol)):>12} {str(np.allclose(x, ref_x, atol=tol)):>12} {str(np.allclose(y, ref_y, atol=tol)):>12} {str(np.allclose(cp, ref_cp, atol=tol)):>12} {str(np.allclose(ue, ref_ue, atol=tol)):>12} {str(np.allclose(dstr, ref_dstar, atol=tol)):>12} {str(np.allclose(theta, ref_theta, atol=tol)):>12} {str(np.allclose(cf, ref_cf, atol=tol)):>12} {str(np.allclose(h, ref_h, atol=tol)):>12} {str(np.allclose(hstar, ref_hstar, atol=tol)):>12}")
        return
 
    elif Description == 'cl':
        a, cl, cd, cm, xtt, xtb, conv = xf.scl(input_start)
        conv_statment = "Converged" if conv else "Not Converged"
        print("")
        print(f"Test {ID}: {Description} - {conv_statment}")
        print("="*115)
        print(f"{'f2pyfoil Scalars:' :<20}" f"A: {a:.10f}, CD: {cd:.10f}, CM: {cm:.10f}, Xtc Top: {xtt:.10f}, Xtc Bot: {xtb:.10f}")
        print(f"{'Reference Scalars:' :<20}" f"A: {ref_alpha[0]:.10f}, CD: {ref_cd[0]:.10f}, CM: {ref_cm[0]:.10f}, Xtc Top: {ref_xtt[0]:.10f}, Xtc Bot: {ref_xtb[0]:.10f}")
        print(f"{'Difference:' :<20}" f"A: {a - ref_alpha[0]:10.5e}, CD: {cd - ref_cd[0]:10.5e}, CM: {cm - ref_cm[0]:10.5e}, Xtc Top: {xtt - ref_xtt[0]:10.5e}, Xtc Bot: {xtb - ref_xtb[0]:10.5e}")
 
        N = xf.getn(); NBL = xf.getnbl()
        cp = xf.getcp(N); cf = xf.getcf(N)
        s, x, y, ue, dstr, theta, tstar, h, hstari = xf.getmorebl(NBL)
        hstar = hstari[: N]
 
        if Re == 0.0:
            print("Cp Data Summary:")
            print(f"{'':>18}{'cp':>12}")
            print(f"{'Max Difference:':<18} {np.max(np.abs(cp - ref_cp)):>12.5e}")
            print(f"{'Mean Difference:':<18} {np.mean(np.abs(cp - ref_cp)):>12.5e}")
            print(f"{'Tolerance Check:':<18} {str(np.allclose(cp, ref_cp, atol=tol)):>12}")
        else:
            print("BL Data Summary:")
            print("="*115)
            print(f"{'':>18} {'s':>12} {'x':>12} {'y':>12} {'Cp':>12} {'Ue/Vinf':>12} {'dstr':>12} {'theta':>12} {'Cf':>12} {'H':>12} {'H*':>12}")
            print(f"{'Max Difference:':<18} {np.max(np.abs(s - ref_s)):>12.5e} {np.max(np.abs(x - ref_x)):>12.5e} {np.max(np.abs(y - ref_y)):>12.5e} {np.max(np.abs(cp - ref_cp)):>12.5e} {np.max(np.abs(ue - ref_ue)):>12.5e} {np.max(np.abs(dstr - ref_dstar)):>12.5e} {np.max(np.abs(theta - ref_theta)):>12.5e} {np.max(np.abs(cf - ref_cf)):>12.5e} {np.max(np.abs(h - ref_h)):>12.5e} {np.max(np.abs(hstar - ref_hstar)):>12.5e}")
            print(f"{'Mean Difference:':<18} {np.mean(np.abs(s - ref_s)):>12.5e} {np.mean(np.abs(x - ref_x)):>12.5e} {np.mean(np.abs(y - ref_y)):>12.5e} {np.mean(np.abs(cp - ref_cp)):>12.5e} {np.mean(np.abs(ue - ref_ue)):>12.5e} {np.mean(np.abs(dstr - ref_dstar)):>12.5e} {np.mean(np.abs(theta - ref_theta)):>12.5e} {np.mean(np.abs(cf - ref_cf)):>12.5e} {np.mean(np.abs(h - ref_h)):>12.5e} {np.mean(np.abs(hstar - ref_hstar)):>12.5e}")
            print(f"{'Tolerance Check:':<18} {str(np.allclose(s, ref_s, atol=tol)):>12} {str(np.allclose(x, ref_x, atol=tol)):>12} {str(np.allclose(y, ref_y, atol=tol)):>12} {str(np.allclose(cp, ref_cp, atol=tol)):>12} {str(np.allclose(ue, ref_ue, atol=tol)):>12} {str(np.allclose(dstr, ref_dstar, atol=tol)):>12} {str(np.allclose(theta, ref_theta, atol=tol)):>12} {str(np.allclose(cf, ref_cf, atol=tol)):>12} {str(np.allclose(h, ref_h, atol=tol)):>12} {str(np.allclose(hstar, ref_hstar, atol=tol)):>12}")
        return
 
    elif Description == 'aseq':
        da = min(np.diff(ref_alpha))
        n_full = int(((input_end - input_start) / da) + 1)
        a, cls, cds, cms, xtts, xtbs, convs = xf.aseq(input_start, input_end, n_full)
        conv_total = sum(convs)
        print("")
        print(f"Test {ID}: {Description} - {conv_total} / {n_full} converged, expected {n_points} / {n_full}")

        mask = np.array([np.any(np.isclose(ai, ref_alpha, atol=tol)) for ai in a])
        print(f"Matching converged simulations: {mask.sum()} / {n_full} matched to reference")
        a    = a[mask];    cls  = cls[mask];  cds  = cds[mask]
        cms  = cms[mask];  xtts = xtts[mask]; xtbs = xtbs[mask]

        print("="*115)
        print(f"{'':>20} {'CL':>10} {'CD':>10} {'CM':>10} {'Xtc Top':>10} {'Xtc Bot':>10}")
        print(f"{'Max Difference:':<20} {np.max(np.abs(cls - ref_cl)):>10.5e} {np.max(np.abs(cds - ref_cd)):>10.5e} {np.max(np.abs(cms - ref_cm)):>10.5e} {np.max(np.abs(xtts - ref_xtt)):>10.5e} {np.max(np.abs(xtbs - ref_xtb)):>10.5e}")
        print(f"{'Mean Difference:':<20} {np.mean(np.abs(cls - ref_cl)):>10.5e} {np.mean(np.abs(cds - ref_cd)):>10.5e} {np.mean(np.abs(cms - ref_cm)):>10.5e} {np.mean(np.abs(xtts - ref_xtt)):>10.5e} {np.mean(np.abs(xtbs - ref_xtb)):>10.5e}")
        print(f"{'Tolerance Check:':<20}", np.allclose(cls, ref_cl, atol=tol), np.allclose(cds, ref_cd, atol=tol), np.allclose(cms, ref_cm, atol=tol), np.allclose(xtts, ref_xtt, atol=tol), np.allclose(xtbs, ref_xtb, atol=tol))
 
    elif Description == 'cseq':
        dcl = min(np.diff(ref_cl))
        n_full = int(((input_end - input_start) / dcl) + 1)
        a, cls, cds, cms, xtts, xtbs, convs = xf.cseq(input_start, input_end + dcl, n_full)
        conv_total = sum(convs)
        print("")
        print(f"Test {ID}: {Description} - {conv_total} / {n_full} converged, expected {n_points} / {n_full}")

        mask = np.array([np.any(np.isclose(cli, ref_cl, atol=tol)) for cli in cls])
        print(f"Matching converged simulations: {mask.sum()} / {n_full} matched to reference")
        a    = a[mask];    cls  = cls[mask];  cds  = cds[mask]
        cms  = cms[mask];  xtts = xtts[mask]; xtbs = xtbs[mask]

        print("="*115)
        print(f"{'':>20} {'Alpha':>10} {'CD':>10} {'CM':>10} {'Xtc Top':>10} {'Xtc Bot':>10}")
        print(f"{'Max Difference:':<20} {np.max(np.abs(a - ref_alpha)):>10.5e} {np.max(np.abs(cds - ref_cd)):>10.5e} {np.max(np.abs(cms - ref_cm)):>10.5e} {np.max(np.abs(xtts - ref_xtt)):>10.5e} {np.max(np.abs(xtbs - ref_xtb)):>10.5e}")
        print(f"{'Mean Difference:':<20} {np.mean(np.abs(a - ref_alpha)):>10.5e} {np.mean(np.abs(cds - ref_cd)):>10.5e} {np.mean(np.abs(cms - ref_cm)):>10.5e} {np.mean(np.abs(xtts - ref_xtt)):>10.5e} {np.mean(np.abs(xtbs - ref_xtb)):>10.5e}")
        print(f"{'Tolerance Check:':<20}", np.allclose(a, ref_alpha, atol=tol), np.allclose(cds, ref_cd, atol=tol), np.allclose(cms, ref_cm, atol=tol), np.allclose(xtts, ref_xtt, atol=tol), np.allclose(xtbs, ref_xtb, atol=tol))
 
Test("test_data.db", 0) 
Test("test_data.db", 1)
Test("test_data.db", 2)
Test("test_data.db", 3)
Test("test_data.db", 4)
Test("test_data.db", 5)
Test("test_data.db", 6)
Test("test_data.db", 7)
Test("test_data.db", 8)
Test("test_data.db", 9)
Test("test_data.db", 10)