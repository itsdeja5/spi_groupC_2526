#!/usr/bin/env python3
"""
best_decaps.py
==============
Finds the decoupling-capacitor configuration (6 values on 5 slots = 7776 cases) that
minimizes the peak-to-peak ripple of the die voltage, starting from LTspice output files.
It writes ONE file, best_config.txt, with the best configuration ready to be pasted
into pdn_best.cir.

INPUT FILES  (in the folder given with --dir, default: current folder)
    zinj_c01 ... zinj_c05, zinj_vdie   1 A injected at one port (no decaps): column of the 6x6 impedance matrix Z(f)
    caps_zgrid                         impedance Z_cap(f) of the 6 capacitors (traces V(z103) ... V(z152))
    die_model/die_model.sp             die current waveform (PWL)
  The extension is .txt (LTspice "Export data as Text", default) or .raw (--format raw). In the .txt files
  every trace must be exported as complex numbers: either (real, imaginary) pairs, as in
  3.959e-002,-2.529e-001   or (magnitude dB, phase deg) pairs, as in  (-1.18e+001dB,-8.11e+001°).
  The column names must be V(c01) ... V(vdie) / V(z103) ... V(z152).
  All files use the same 59999 frequencies  f_k = k * df,  df = 1/(N*dt)  (FFT of the die current).

THE MATH  (all quantities are complex, one value per frequency f_k)
  Z is the 6x6 impedance matrix of the PDN. Ports: 5 capacitor ports (c) and the die port (d):

      V_c = Z_cc * I_c + Z_cd * I_d        (1)   currents I flow INTO the network
      V_d = Z_dc * I_c + Z_dd * I_d        (2)

  A capacitor connected to a port takes the current that leaves the network, and obeys
  its own Ohm's law (Z_cap = diagonal matrix with the 5 chosen capacitor impedances):

      V_c = - Z_cap * I_c                  (3)

  Replacing (3) in (1) gives I_c = -(Z_cc + Z_cap)^-1 * Z_cd * I_d, and then (2) gives
  the impedance seen by the die with the capacitors connected:

      Z_die(f) = Z_dd - Z_dc * (Z_cc + Z_cap)^-1 * Z_cd                       <-- (A)

  The die draws the known current I_load(f) (FFT of the PWL waveform), so the ripple
  of the die voltage, at every frequency, is

      V_die(f) = - Z_die(f) * I_load(f)                                       <-- (B)

  (minus sign: the die DRAWS current from the supply). The inverse FFT of (B) is the
  ripple waveform v_die(t); the figure of merit is  Vpp = max(v_die) - min(v_die).
  The DC drop is ignored: it does not change the peak-to-peak value.
  The script evaluates (A) and (B) for the 7776 configurations and keeps the minimum Vpp.

USAGE
    pip install numpy            (spicelib too, only if you use --format raw)
    python best_decaps.py                    # files in the current folder
    python best_decaps.py --dir SW2_STM      # files in another folder
    python best_decaps.py --format raw       # if you use the .raw files
"""
import argparse
import itertools
import os
import re
import sys
import time

import numpy as np

# ----------------------------------------------------------------------------------
# FILE NAMES AND DATA  (edit here if your names are different)
# ----------------------------------------------------------------------------------
PORTS = ['c01', 'c02', 'c03', 'c04', 'c05', 'vdie']     # order of rows/columns of Z(f)
MATRIX_BASENAMES = {p: f'zinj_{p}' for p in PORTS}        # 1 A injected at port p (extension added later)
CAPS_BASENAME = 'caps_zgrid'
DIE_FILE = os.path.join('../input_files/die_model', 'die_model.sp')
# capacitors: (code in the part number, label, nominal capacitance [F]); node in caps_zgrid.cir = 'z' + code
CAPS = [('103', '10nF', 10e-9), ('682', '6.8nF', 6.8e-9), ('472', '4.7nF', 4.7e-9),
        ('332', '3.3nF', 3.3e-9), ('222', '2.2nF', 2.2e-9), ('152', '1.5nF', 1.5e-9)]
SLOTS = ['C01', 'C02', 'C03', 'C04', 'C05']
part = lambda code: f'GCM033R71A{code}KA03'
N_SLOTS = 5


# ----------------------------------------------------------------------------------
# READING THE FILES
# ----------------------------------------------------------------------------------
def read_raw_voltages(path, nodes, dialect=None):
    """Frequencies [Hz] and complex node voltages (K x len(nodes)) from an LTspice .raw AC file."""
    try:
        from spicelib import RawRead
    except ImportError:
        sys.exit('The .raw files are read with the spicelib package:  pip install spicelib')
    if not os.path.exists(path):
        sys.exit(f'file not found: {path}')
    try:
        raw = RawRead(path, dialect=dialect) if dialect else RawRead(path)
    except Exception:
        if dialect:
            raise
        raw = RawRead(path, dialect='ngspice')       # only needed for files not written by LTspice
    names = {n.lower(): n for n in raw.get_trace_names()}
    try:                                                  # LTspice: frequency is the axis of the file
        f = np.real(raw.get_axis())
    except Exception:                                     # other simulators: frequency is a normal trace
        f = np.real(raw.get_trace(names['frequency']).get_wave())
    cols = []
    for node in nodes:
        key = f'v({node})' if f'v({node})' in names else node
        if key not in names:
            sys.exit(f'{path}: trace V({node}) not found; available: {sorted(names.values())}')
        w = np.asarray(raw.get_trace(names[key]).get_wave())
        if not np.iscomplexobj(w):
            sys.exit(f'{path}: V({node}) is not complex: was it produced by an .ac analysis?')
        cols.append(w.astype(complex))
    return np.abs(f), np.array(cols).T


def read_txt_voltages(path, nodes):
    """Frequencies [Hz] and complex node voltages (K x len(nodes)) from an LTspice 'Export data as Text' file."""
    if not os.path.exists(path):
        sys.exit(f'file not found: {path}')
    data = open(path, 'rb').read()
    if data[:2] in (b'\xff\xfe', b'\xfe\xff'):
        text = data.decode('utf-16')
    else:
        try:
            text = data.decode('utf-8')
        except UnicodeDecodeError:
            text = data.decode('cp1252', errors='replace')
    lines = [l for l in text.splitlines() if l.strip()]
    names = [h.strip() for h in lines[0].split('\t')]            # first column = frequency, then the traces
    body = lines[1:]
    f = np.array([float(l.split('\t', 1)[0]) for l in body])
    joined = '\n'.join(l.split('\t', 1)[1] for l in body)
    if 'dB' in body[0]:                                          # (magnitude dB, phase deg) pairs
        tok = re.findall(r'\(\s*([-+0-9.eE]+)\s*dB\s*,\s*([-+0-9.eE]+)[^)]*\)', joined)
        arr = np.array(tok, dtype=float).reshape(len(body), -1, 2) if tok else None
        C = None if arr is None else 10.0 ** (arr[:, :, 0] / 20.0) * np.exp(1j * np.radians(arr[:, :, 1]))
    else:                                                        # (real, imaginary) pairs  "re,im"  (with or without parentheses)
        num = r'([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)'
        tok = re.findall(num + r'\s*,\s*' + num, joined)
        arr = np.array(tok, dtype=float).reshape(len(body), -1, 2) if tok else None
        C = None if arr is None else arr[:, :, 0] + 1j * arr[:, :, 1]
    if C is None or C.shape[1] != len(names) - 1:
        sys.exit(f'{path}: cannot read the traces as complex pairs. In the export dialog choose the complex '
                 f'format ((dB, phase) or (real, imaginary)); header: {names}')
    low = [n.lower() for n in names]
    cols = []
    for node in nodes:
        key = f'v({node})'
        if key not in low:
            sys.exit(f'{path}: trace {key.upper()} not found; the columns are {names}')
        cols.append(C[:, low.index(key) - 1])
    return f, np.array(cols).T


def read_voltages(path, nodes, fmt, dialect=None):
    return read_txt_voltages(path, nodes) if fmt == 'txt' else read_raw_voltages(path, nodes, dialect)


def read_die_current(path):
    """Die current i(t) [A] from the PWL source in die_model.sp (times in ps, uniform step)."""
    txt = open(path).read().replace('\r', '')
    body = txt.split('pwl(')[1].split(')')[0]
    pairs = re.findall(r'([0-9.eE+-]+)ps\s+([0-9.eE+-]+)', body)
    t = np.array([float(a) for a, _ in pairs]) * 1e-12
    i = np.array([float(b) for _, b in pairs])
    dt = t[1] - t[0]
    if not np.allclose(np.diff(t), dt, rtol=1e-6):
        sys.exit('the PWL current is not uniformly sampled')
    return i, dt


# ----------------------------------------------------------------------------------
# THE TWO KEY EXPRESSIONS:  (A) Z_die  and  (B) V_die
# ----------------------------------------------------------------------------------
def die_impedance(Z, Z_cap_values, lines):
    """
    (A)   Z_die(f) = Z_dd - Z_dc * (Z_cc + Z_cap)^-1 * Z_cd
    Z            : (K, 6, 6) impedance matrix of the PDN without capacitors, ports ordered C01..C05, die
    Z_cap_values : (len(lines), 5) impedance of the capacitor mounted on each of the 5 slots
    lines        : indices of the frequencies to evaluate
    """
    Z_cc = Z[lines, :N_SLOTS, :N_SLOTS].copy()          # 5x5  capacitor ports <-> capacitor ports
    Z_cd = Z[lines, :N_SLOTS, N_SLOTS]                  # 5x1  capacitor ports <- die port
    Z_dc = Z[lines, N_SLOTS, :N_SLOTS]                  # 1x5  die port <- capacitor ports
    Z_dd = Z[lines, N_SLOTS, N_SLOTS]                   # 1x1  die port <- die port (no caps)
    idx = np.arange(N_SLOTS)
    Z_cc[:, idx, idx] += Z_cap_values                   # Z_cc + Z_cap   (Z_cap is diagonal)
    x = np.linalg.solve(Z_cc, Z_cd[:, :, None])[:, :, 0]        # (Z_cc + Z_cap)^-1 * Z_cd
    return Z_dd - np.einsum('ki,ki->k', Z_dc, x)                # Z_dd - Z_dc * x


def die_ripple_waveform(Z_die, I_load, lines, n):
    """
    (B)   V_die(f) = - Z_die(f) * I_load(f)       then inverse FFT  ->  v_die(t)
    I_load : FFT of the die current, line k at f_k = k*df  (index k-1 in the arrays)
    n      : number of time samples of the die current
    """
    V_die = np.zeros(n // 2 + 1, complex)               # spectrum of the ripple (FFT bins 0..n/2)
    V_die[lines + 1] = -Z_die * I_load[lines]           # <-- V_die(f_k) = -Z_die(f_k) * I_load(f_k)
    return np.fft.irfft(V_die, n)                       # ripple waveform v_die(t) [V]


def ripple_vpp(Z, Zc, cfg, I_load, lines, n):
    """Peak-to-peak ripple [V] of the configuration cfg (5 indices into CAPS, one per slot)."""
    Z_cap_values = Zc[lines][:, list(cfg)]
    v = die_ripple_waveform(die_impedance(Z, Z_cap_values, lines), I_load, lines, n)
    return v.max() - v.min()


# ----------------------------------------------------------------------------------
# CONSISTENCY CHECKS (stop with a clear message if something is wrong with the files)
# ----------------------------------------------------------------------------------
def check_inputs(f, Z, Zc, df):
    K = len(f)
    f_exp = (np.arange(K) + 1) * df
    err = np.max(np.abs(f - f_exp) / f_exp)
    if err > 1e-5:
        sys.exit(f'the frequencies in the files are not k*df (df = {df:.4f} Hz): max relative error {err:.1e}.\n'
                 f'Check the .ac line of the netlists (.ac lin {K} {df!r} {K * df!r}).')
    asym = np.max(np.abs(Z - Z.transpose(0, 2, 1)) / np.abs(Z))
    if asym > 1e-2:
        sys.exit(f'Z(f) is not symmetric (max relative asymmetry {asym:.1e}): the zinj_* files are probably '
                 f'mixed up, or the injection node was not changed correctly in one of the copies.')
    if asym > 1e-5:
        print(f'  warning: Z(f) symmetric only to {asym:.1e}')
    for j, (code, label, C) in enumerate(CAPS):             # low-frequency capacitor impedance ~ 1/(2 pi f C)
        z0 = abs(Zc[0, j])
        z_ideal = 1.0 / (2 * np.pi * f[0] * C)
        if abs(z0 - z_ideal) / z_ideal > 0.25:
            print(f'  warning: |Z_cap| of {label} at {f[0] / 1e3:.1f} kHz is {z0:.1f} ohm, expected about '
                  f'{z_ideal:.1f}: are the traces of the capacitor file in the right order?')
    print(f'  checks passed: grid k*df ({K} lines up to {f[-1] / 1e9:.3f} GHz), Z symmetric to {asym:.1e}')


# ----------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description='best decoupling-capacitor configuration from LTspice outputs')
    ap.add_argument('--dir', default='.', help='folder with the .raw files and die_model/ (default: .)')
    ap.add_argument('--out', default='best_config.txt', help='output file (default: best_config.txt)')
    ap.add_argument('--top', type=int, default=10, help='how many runner-up configurations to list')
    ap.add_argument('--cover', type=float, default=0.999, help='fraction of the die-current power used in the screening')
    ap.add_argument('--nref', type=int, default=300, help='best candidates re-evaluated with all frequencies')
    ap.add_argument('--format', default='txt', choices=['txt', 'raw'], help='type of the LTspice output files (default: txt)')
    ap.add_argument('--dialect', default=None, help='raw-file dialect (leave empty for LTspice files)')
    args = ap.parse_args()
    P = lambda name: os.path.join(args.dir, name)

    # 1) die current and its spectrum I_load(f_k)
    cur, dt = read_die_current(P(DIE_FILE))
    n = len(cur)
    df = 1.0 / (n * dt)
    I_all = np.fft.rfft(cur)
    print(f'die current: {n} samples, dt = {dt * 1e12:.0f} ps, df = {df:.2f} Hz')

    # 2) impedance matrix Z(f): column j = voltages at the 6 ports with 1 A injected at port j
    print('reading the impedance-matrix files ...')
    Z = None
    for j, p in enumerate(PORTS):
        fname = f'{MATRIX_BASENAMES[p]}.{args.format}'
        f, V = read_voltages(P(fname), PORTS, args.format, args.dialect)
        if Z is None:
            f0, Z = f, np.zeros((len(f), 6, 6), complex)
        elif len(f) != len(f0) or np.max(np.abs(f - f0) / f0) > 1e-9:
            sys.exit(f'{fname}: frequency grid different from the other files')
        Z[:, :, j] = V                                    # Z(i, j) = V_i with 1 A at port j
    K = len(f0)

    # 3) capacitor impedances Z_cap(f)
    print('reading the capacitor file ...')
    caps_name = f'{CAPS_BASENAME}.{args.format}'
    fc, Vc = read_voltages(P(caps_name), [f'z{c}' for c, _, _ in CAPS], args.format, args.dialect)
    if len(fc) != K or np.max(np.abs(fc - f0) / f0) > 1e-9:
        sys.exit(f'{caps_name}: frequency grid different from the matrix files')
    Zc = Vc                                               # (K, 6): Z_cap of each of the 6 capacitors
    check_inputs(f0, Z, Zc, df)

    I_load = I_all[1:K + 1]                               # die-current spectrum on the same frequencies
    E = np.abs(I_load) ** 2
    order = np.argsort(E)[::-1]
    ncut = int(np.searchsorted(np.cumsum(E[order]) / E.sum(), args.cover)) + 1
    sel = np.sort(order[:ncut])                           # lines holding `cover` of the current power
    allk = np.arange(K)
    print(f'screening on {ncut} of {K} frequency lines ({args.cover * 100:.1f} % of the current power)')

    # 4) ripple of every configuration; no-decap reference
    vpp_nocaps = die_ripple_waveform(Z[:, N_SLOTS, N_SLOTS], I_load, allk, n)
    vpp_nocaps = vpp_nocaps.max() - vpp_nocaps.min()
    cfgs = list(itertools.product(range(len(CAPS)), repeat=N_SLOTS))
    pk = np.zeros(len(cfgs))
    t0 = time.time()
    for m, c in enumerate(cfgs):
        pk[m] = ripple_vpp(Z, Zc, c, I_load, sel, n)
        if m % 1000 == 0:
            print(f'  {m}/{len(cfgs)} configurations evaluated ({time.time() - t0:.0f} s)')
    top = np.argsort(pk)[:args.nref]
    final = {int(m): ripple_vpp(Z, Zc, cfgs[m], I_load, allk, n) for m in top}      # refine with all lines
    rank = sorted(final, key=final.get)
    print(f'search finished in {time.time() - t0:.0f} s')

    # 5) the single output file
    best = cfgs[rank[0]]
    with open(P(args.out) if args.dir != '.' else args.out, 'w') as fo:
        fo.write('BEST DECOUPLING-CAPACITOR CONFIGURATION (minimum peak-to-peak ripple at the die)\n')
        fo.write('=' * 80 + '\n')
        fo.write('Method:  Z_die(f) = Z_dd - Z_dc * inv(Z_cc + Z_cap) * Z_cd\n')
        fo.write('         V_die(f) = - Z_die(f) * I_load(f)   ->  inverse FFT  ->  Vpp = max - min\n\n')
        fo.write(f'Peak-to-peak ripple without decaps : {vpp_nocaps:.4f} V\n')
        fo.write(f'Peak-to-peak ripple, best          : {final[rank[0]]:.4f} V\n\n')
        fo.write('Best configuration. Paste these lines into pdn_best.cir:\n')
        for s, k in enumerate(best):
            fo.write(f'X{SLOTS[s].lower()} {SLOTS[s].lower()} 0 {part(CAPS[k][0])}      * {CAPS[k][1]}\n')
        fo.write(f'\nNext {args.top} configurations (Vpp [V] | C01  C02  C03  C04  C05):\n')
        for r in rank[1:args.top + 1]:
            fo.write(f'  {final[r]:.4f} | ' + '  '.join(f'{CAPS[k][1]:>5s}' for k in cfgs[r]) + '\n')
        fo.write('\nNotes:\n')
        fo.write('  - Vpp is computed on the periodic steady-state waveform of the 12 us die current.\n')
        fo.write('  - Confirm the chosen configuration with a transient in LTspice (pdn_best.cir) and measure\n')
        fo.write('    the ripple from 0.5 us on, to skip the start-up transient.\n')
    print(f'\nbest: ' + ', '.join(f'{SLOTS[s]}={CAPS[k][1]}' for s, k in enumerate(best))
          + f'  ->  Vpp = {final[rank[0]]:.4f} V   (no decaps: {vpp_nocaps:.4f} V)')
    print(f'written: {args.out}')


if __name__ == '__main__':
    main()
