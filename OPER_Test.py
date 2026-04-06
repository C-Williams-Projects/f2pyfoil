import numpy as np
from OPER_data import (Test_5, Test_6, Test_7, Test_Descriptions, Test_0, Test_1, Test_2, Test_3, Test_4,
                    Test_5, Test_6, Test_7, Test_8, Test_9, DAE_11, PANELS_DAE_11, PANELS_NACA0012, PANELS_NACA2412)
import xfoil as xf

def Single(case, tol=1e-6):
    ID = case['Test_ID']
    Description = Test_Descriptions[str(ID)] 
    xf.initialise()
    xf.quiet(True)
    if case['Airfoil'] == 'DAE-11':
        xf.airfoil(DAE_11['x'], DAE_11['y'])
        ref_x = PANELS_DAE_11['x']
        ref_y = PANELS_DAE_11['y']
        ref_s = PANELS_DAE_11['s']
    else:
        xf.naca_load(case['Airfoil'])
        if case['Airfoil'] == '0012':
            ref_x = PANELS_NACA0012['x']
            ref_y = PANELS_NACA0012['y']
            ref_s = PANELS_NACA0012['s']
        else:
            ref_x = PANELS_NACA2412['x']
            ref_y = PANELS_NACA2412['y']
            ref_s = PANELS_NACA2412['s']
    xf.repanel(160, 1.0, 0.15, 0.2, 1.0, 1.0, 1.0, 1.0)
    xf.setcon(case['Re'], case['Ma'])
    xf.setiter(300)
    xf.trpars(case['XTR Top'], case['XTR Bot'], case['Ncrit'])
    if Description == 'Single alpha':
        cl, cd, cm, xtt, xtb, conv = xf.alpha(case['alpha'])  
        conv_statment = "Converged" if conv else "Not Converged"
        print("")
        print(f"Test Case {ID}: {Description} - {conv_statment}")
        print("="*20)  
        print("f2pyfoil Scalars{:<7}:" f"CL: {cl:.10f}, CD: {cd:.10f}, CM: {cm:.10f}, Xtc Top: {xtt:.10f}, Xtc Bot: {xtb:.10f}")
        print("XFoil Reference Scalars:" f"CL: {case['CL']:.10f}, CD: {case['CD']:.10f}, CM: {case['CM']:.10f}, Xtc Top: {case['Xtc Top']:.10f}, Xtc Bot: {case['Xtc Bot']:.10f}")
        print("Difference in Scalars:" f"CL: {cl - case['CL']:.10f}, CD: {cd - case['CD']:.10f}, CM: {cm - case['CM']:.10f}, Xtc Top: {xtt - case['Xtc Top']:.10f}, Xtc Bot: {xtb - case['Xtc Bot']:.10f}")
    elif Description == 'Single cl':
        a, cd, cm, xtt, xtb, conv = xf.cl(case['cl'])
        conv_statment = "Converged" if conv else "Not Converged"
        print("")
        print(f"Test Case {ID}: {Description} - {conv_statment}") 
        print("="*20) 
        print("f2pyfoil Scalars{:<7}:" f"A: {a:.10f}, CD: {cd:.10f}, CM: {cm:.10f}, Xtc Top: {xtt:.10f}, Xtc Bot: {xtb:.10f}")
        print("XFoil Reference Scalars:" f"A: {case['CL']:.10f}, CD: {case['CD']:.10f}, CM: {case['CM']:.10f}, Xtc Top: {case['Xtc Top']:.10f}, Xtc Bot: {case['Xtc Bot']:.10f}")
        print("Difference in Scalars:" f"A: {a - case['alpha']:.10f}, CD: {cd - case['CD']:.10f}, CM: {cm - case['CM']:.10f}, Xtc Top: {xtt - case['Xtc Top']:.10f}, Xtc Bot: {xtb - case['Xtc Bot']:.10f}")
    else:
        cl, cd, cm, xtt, xtb, conv = xf.alpha(0.0)
        cl, cd, cm, xtt, xtb, conv = xf.alpha(case['alpha'] / 4)
        cl, cd, cm, xtt, xtb, conv = xf.alpha(case['alpha'] / 3)
        cl, cd, cm, xtt, xtb, conv = xf.alpha(case['alpha'] / 2)
        cl, cd, cm, xtt, xtb, conv = xf.alpha(2 * case['alpha'] / 3)
        cl, cd, cm, xtt, xtb, conv = xf.alpha(case['alpha'])
        conv_statment = "Converged" if conv else "Not Converged"

        print("")
        print(f"Test Case {ID}: {Description} - {conv_statment}") 
        print("="*20)  
        print("f2pyfoil Scalars{:<7}:" f"CL: {cl:.10f}, CD: {cd:.10f}, CM: {cm:.10f}, Xtc Top: {xtt:.10f}, Xtc Bot: {xtb:.10f}")
        print("XFoil Reference Scalars:" f"CL: {case['CL']:.10f}, CD: {case['CD']:.10f}, CM: {case['CM']:.10f}, Xtc Top: {case['Xtc Top']:.10f}, Xtc Bot: {case['Xtc Bot']:.10f}")
        print("Difference in Scalars:" f"CL: {cl - case['CL']:.10f}, CD: {cd - case['CD']:.10f}, CM: {cm - case['CM']:.10f}, Xtc Top: {xtt - case['Xtc Top']:.10f}, Xtc Bot: {xtb - case['Xtc Bot']:.10f}")
        
    N = xf.getn(); NBL = xf.getnbl()
    cp = xf.getcp(N)
    s, xi, y, ue, dstr, thet, cf, hk, cdis, ctau, hstar, P, m, K  = xf.getbl(NBL)
    x = xi[:N]
    if case['Re'] == 0.0:
        print("Cp Data Summary:")
        print("="*45)
        print(f"{'':>18} {'x':>12} {'cp':>12}")
        print(f"{'Max Difference:':>18} {np.max(np.abs(x - ref_x)):>12.5e} {np.max(np.abs(cp - case['Cp'])):>12.5e}")
        print(f"{'Mean Difference:':>18} {np.mean(np.abs(x - ref_x)):>12.5e} {np.mean(np.abs(cp - case['Cp'])):>12.5e}")
        print(f"{'Tolerance Check:':>18} {str(np.allclose(x, ref_x, atol=tol)):>12} {str(np.allclose(cp, case['Cp'], atol=tol)):>12}")
    else:
        print("BL Data Summary:")
        print("="*115)
        print(f"{'':>18} {'s':>12} {'x':>12} {'y':>12} {'cp':>12} {'ue':>12} {'dstr':>12} {'thet':>12} {'cf':>12}")
        print(f"{'Max Difference:':>18} {np.max(np.abs(s - ref_s)):>12.5e} {np.max(np.abs(x - ref_x)):>12.5e} {np.max(np.abs(y - ref_y)):>12.5e} {np.max(np.abs(cp - case['Cp'])):>12.5e} {np.max(np.abs(ue - case['Ue'])):>12.5e} {np.max(np.abs(dstr - case['Dstar'])):>12.5e} {np.max(np.abs(thet - case['Theta'])):>12.5e} {np.max(np.abs(cf - case['Cf'])):>12.5e}")
        print(f"{'Mean Difference:':>18} {np.mean(np.abs(s - ref_s)):>12.5e} {np.mean(np.abs(x - ref_x)):>12.5e} {np.mean(np.abs(y - ref_y)):>12.5e} {np.mean(np.abs(cp - case['Cp'])):>12.5e} {np.mean(np.abs(ue - case['Ue'])):>12.5e} {np.mean(np.abs(dstr - case['Dstar'])):>12.5e} {np.mean(np.abs(thet - case['Theta'])):>12.5e} {np.mean(np.abs(cf - case['Cf'])):>12.5e}")
        print(f"{'Tolerance Check:':>18} {str(np.allclose(s, ref_s, atol=tol)):>12} {str(np.allclose(x, ref_x, atol=tol)):>12} {str(np.allclose(y, ref_y, atol=tol)):>12} {str(np.allclose(cp, case['Cp'], atol=tol)):>12} {str(np.allclose(ue, case['Ue'], atol=tol)):>12} {str(np.allclose(dstr, case['Dstar'], atol=tol)):>12} {str(np.allclose(thet, case['Theta'], atol=tol)):>12} {str(np.allclose(cf, case['Cf'], atol=tol)):>12}")
    return

def SingleWarmup(case, tol=1e-6):
    ID = case['Test_ID']
    Description = Test_Descriptions[str(ID)] 
    xf.initialise()
    xf.quiet(True)
    if case['Airfoil'] == 'DAE-11':
        xf.airfoil(DAE_11['x'], DAE_11['y'])
    else:
        xf.naca_load(case['Airfoil'])
    xf.repanel(160, 1.0, 0.15, 0.2, 1.0, 1.0, 1.0, 1.0)
    xf.setcon(case['Re'], case['Ma'])
    xf.setiter(300)
    xf.trpars(case['XTR Top'], case['XTR Bot'], case['Ncrit'])
    if Description == 'Single alpha':
        cl, cd, cm, xtt, xtb, conv = xf.alpha(case['alpha'])  
        conv_statment = "Converged" if conv else "Not Converged"
        print("")
        print(f"Warmup Test Case {ID}: {Description} - {conv_statment}")
        print("="*20)  
        print("f2pyfoil Scalars{:<7}:" f"CL: {cl:.10f}, CD: {cd:.10f}, CM: {cm:.10f}, Xtc Top: {xtt:.10f}, Xtc Bot: {xtb:.10f}")
        print("XFoil Reference Scalars:" f"CL: {case['CL']:.10f}, CD: {case['CD']:.10f}, CM: {case['CM']:.10f}, Xtc Top: {case['Xtc Top']:.10f}, Xtc Bot: {case['Xtc Bot']:.10f}")


def Sequential(case, tol=1e-6):
    ID = case['Test_ID']
    Description = Test_Descriptions[str(ID)] 
    xf.initialise()
    xf.quiet(True)
    if case['Airfoil'] == 'DAE-11':
        xf.airfoil(DAE_11['x'], DAE_11['y'])
    else:
        xf.naca_load(case['Airfoil'])
    xf.repanel(160, 1.0, 0.15, 0.2, 1.0, 1.0, 1.0, 1.0)
    xf.setcon(case['Re'], case['Ma'])
    xf.setiter(300)
    xf.trpars(case['XTR Top'], case['XTR Bot'], case['Ncrit'])
    if Description == 'Single alpha':
        cls, cds, cms, xtts, xtbs, convs = xf.aseq(case['alpha_i'], case['alpha_f'], case['alpha_step'])
        total_converged = sum(convs)
        print("")
        print(f"Test Case {ID}: {Description} - {total_converged}/{len(convs)} Converged")
        print("="*20)
        print(f"{'CL':>10} {'CD':>10} {'CM':>10} {'Xtc Top':>10} {'Xtc Bot':>10}")
        print("Max Difference:", f"{np.max(np.abs(cls - case['CLs'])):>10.5e} {np.max(np.abs(cds - case['CDs'])):>10.5e} {np.max(np.abs(cms - case['CMs'])):>10.5e} {np.max(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.max(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Mean Difference:", f"{np.mean(np.abs(cls - case['CLs'])):>10.5e} {np.mean(np.abs(cds - case['CDs'])):>10.5e} {np.mean(np.abs(cms - case['CMs'])):>10.5e} {np.mean(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.mean(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Tolerance Check:", np.allclose(cls, case['CLs'], atol=tol), np.allclose(cds, case['CDs'], atol=tol), np.allclose(cms, case['CMs'], atol=tol), np.allclose(xtts, case['Xtc Tops'], atol=tol), np.allclose(xtbs, case['Xtc Bots'], atol=tol))

    else:
        alphas, cds, cms, xtts, xtbs, convs = xf.cseq(case['cl_i'], case['cl_f'], case['cl_step'])
        total_converged = sum(convs)
        print("")
        print(f"Test Case {ID}: {Description} - {total_converged}/{len(convs)} Converged")
        print("="*20)
        print(f"{'CL':>10} {'CD':>10} {'CM':>10} {'Xtc Top':>10} {'Xtc Bot':>10}")
        print("Max Difference:", f"{np.max(np.abs(alphas - case['alphas'])):>10.5e} {np.max(np.abs(cds - case['CDs'])):>10.5e} {np.max(np.abs(cms - case['CMs'])):>10.5e} {np.max(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.max(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Mean Difference:", f"{np.mean(np.abs(alphas - case['alphas'])):>10.5e} {np.mean(np.abs(cds - case['CDs'])):>10.5e} {np.mean(np.abs(cms - case['CMs'])):>10.5e} {np.mean(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.mean(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Tolerance Check:", np.allclose(alphas, case['alphas'], atol=tol), np.allclose(cds, case['CDs'], atol=tol), np.allclose(cms, case['CMs'], atol=tol), np.allclose(xtts, case['Xtc Tops'], atol=tol), np.allclose(xtbs, case['Xtc Bots'], atol=tol))
   

def Consecutive(case, tol=1e-6):
    ID = case['Test_ID']
    Description = Test_Descriptions[str(ID)] 
    xf.initialise()
    xf.quiet(True)
    if case['Airfoil'] == 'DAE-11':
        xf.airfoil(DAE_11['x'], DAE_11['y'])
        ref_x = PANELS_DAE_11['x']
        ref_y = PANELS_DAE_11['y']
        ref_s = PANELS_DAE_11['s']
    else:
        xf.naca_load(case['Airfoil'])
        if case['Airfoil'] == '0012':
            ref_x = PANELS_NACA0012['x']
            ref_y = PANELS_NACA0012['y']
            ref_s = PANELS_NACA0012['s']
        else:
            ref_x = PANELS_NACA2412['x']
            ref_y = PANELS_NACA2412['y']
            ref_s = PANELS_NACA2412['s']
    xf.repanel(160, 1.0, 0.15, 0.2, 1.0, 1.0, 1.0, 1.0)
    xf.setcon(case['Re'], case['Ma'])
    xf.setiter(300)
    xf.trpars(case['XTR Top'], case['XTR Bot'], case['Ncrit'])
    if Description == 'Consecutive alpha':
        cls, cds, cms, xtts, xtbs, convs = [], [], [], [], [], []
        for i in range(len(case['alphas'])):
            cl, cd, cm, xtt, xtb, conv = xf.alpha(case['alphas'][i])
            cls.append(cl); cds.append(cd); cms.append(cm)
            xtts.append(xtt); xtbs.append(xtb); convs.append(conv)    
        total_converged = sum(convs)
        print(f"Test Case {ID}: {Description} - {total_converged}/{len(convs)} Converged")
        print(f"{'CL':>10} {'CD':>10} {'CM':>10} {'Xtc Top':>10} {'Xtc Bot':>10}")
        print("Max Difference:", f"{np.max(np.abs(cls - case['CLs'])):>10.5e} {np.max(np.abs(cds - case['CDs'])):>10.5e} {np.max(np.abs(cms - case['CMs'])):>10.5e} {np.max(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.max(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Mean Difference:", f"{np.mean(np.abs(cls - case['CLs'])):>10.5e} {np.mean(np.abs(cds - case['CDs'])):>10.5e} {np.mean(np.abs(cms - case['CMs'])):>10.5e} {np.mean(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.mean(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Tolerance Check:", np.allclose(cls, case['CLs'], atol=tol), np.allclose(cds, case['CDs'], atol=tol), np.allclose(cms, case['CMs'], atol=tol), np.allclose(xtts, case['Xtc Tops'], atol=tol), np.allclose(xtbs, case['Xtc Bots'], atol=tol))
        print("")
    else:
        alphas, cds, cms, xtts, xtbs, convs = [], [], [], [], [], []
        for i in range(len(case['cls'])):
            alpha, cd, cm, xtt, xtb, conv = xf.cl(case['cls'][i])
            alphas.append(alpha); cds.append(cd); cms.append(cm)
            xtts.append(xtt); xtbs.append(xtb); convs.append(conv)
        total_converged = sum(convs)
        print("")
        print(f"Test Case {ID}: {Description} - {total_converged}/{len(convs)} Converged")
        print("="*20)
        print(f"{'Alpha':>10} {'CD':>10} {'CM':>10} {'Xtc Top':>10} {'Xtc Bot':>10}")
        print("Max Difference:", f"{np.max(np.abs(alphas - case['alphas'])):>10.5e} {np.max(np.abs(cds - case['CDs'])):>10.5e} {np.max(np.abs(cms - case['CMs'])):>10.5e} {np.max(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.max(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Mean Difference:", f"{np.mean(np.abs(alphas - case['alphas'])):>10.5e} {np.mean(np.abs(cds - case['CDs'])):>10.5e} {np.mean(np.abs(cms - case['CMs'])):>10.5e} {np.mean(np.abs(xtts - case['Xtc Tops'])):>10.5e} {np.mean(np.abs(xtbs - case['Xtc Bots'])):>10.5e}")
        print("Tolerance Check:", np.allclose(alphas, case['alphas'], atol=tol), np.allclose(cds, case['CDs'], atol=tol), np.allclose(cms, case['CMs'], atol=tol), np.allclose(xtts, case['Xtc Tops'], atol=tol), np.allclose(xtbs, case['Xtc Bots'], atol=tol))
        print("")
 
    N = xf.getn(); NBL = xf.getnbl()
    cp = xf.getcp(N)
    s, xi, y, ue, dstr, thet, cf, hk, cdis, ctau, hstar, P, m, K  = xf.getbl(NBL)
    x = xi[:N]
    if case['Re'] == 0.0:
        print("Cp Data Summary:")
        print("="*45)
        print(f"{'':>18} {'x':>12} {'cp':>12}")
        print(f"{'Max Difference:':>18} {np.max(np.abs(x - ref_x)):>12.5e} {np.max(np.abs(cp - case['Cp'])):>12.5e}")
        print(f"{'Mean Difference:':>18} {np.mean(np.abs(x - ref_x)):>12.5e} {np.mean(np.abs(cp - case['Cp'])):>12.5e}")
        print(f"{'Tolerance Check:':>18} {str(np.allclose(x, ref_x, atol=tol)):>12} {str(np.allclose(cp, case['Cp'], atol=tol)):>12}")
    else:
        print("BL Data Summary:")
        print("="*115)
        print(f"{'':>18} {'s':>12} {'x':>12} {'y':>12} {'cp':>12} {'ue':>12} {'dstr':>12} {'thet':>12} {'cf':>12}")
        print(f"{'Max Difference:':>18} {np.max(np.abs(s - ref_s)):>12.5e} {np.max(np.abs(x - ref_x)):>12.5e} {np.max(np.abs(y - ref_y)):>12.5e} {np.max(np.abs(cp - case['Cp'])):>12.5e} {np.max(np.abs(ue - case['Ue'])):>12.5e} {np.max(np.abs(dstr - case['Dstar'])):>12.5e} {np.max(np.abs(thet - case['Theta'])):>12.5e} {np.max(np.abs(cf - case['Cf'])):>12.5e}")
        print(f"{'Mean Difference:':>18} {np.mean(np.abs(s - ref_s)):>12.5e} {np.mean(np.abs(x - ref_x)):>12.5e} {np.mean(np.abs(y - ref_y)):>12.5e} {np.mean(np.abs(cp - case['Cp'])):>12.5e} {np.mean(np.abs(ue - case['Ue'])):>12.5e} {np.mean(np.abs(dstr - case['Dstar'])):>12.5e} {np.mean(np.abs(thet - case['Theta'])):>12.5e} {np.mean(np.abs(cf - case['Cf'])):>12.5e}")
        print(f"{'Tolerance Check:':>18} {str(np.allclose(s, ref_s, atol=tol)):>12} {str(np.allclose(x, ref_x, atol=tol)):>12} {str(np.allclose(y, ref_y, atol=tol)):>12} {str(np.allclose(cp, case['Cp'], atol=tol)):>12} {str(np.allclose(ue, case['Ue'], atol=tol)):>12} {str(np.allclose(dstr, case['Dstar'], atol=tol)):>12} {str(np.allclose(thet, case['Theta'], atol=tol)):>12} {str(np.allclose(cf, case['Cf'], atol=tol)):>12}")

Single(Test_0)
Single(Test_1)
Single(Test_2)
Single(Test_3)
Single(Test_4)