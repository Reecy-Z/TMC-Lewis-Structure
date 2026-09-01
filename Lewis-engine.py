"""
NBO-style Best Lewis Structure Finder
======================================
Finds the best Lewis structure and prints the $CHOOSE block,
replicating the exact NBO algorithm from nbo.log.

Algorithm
---------
1. CONNECTIVITY
   Bond if dist(i,j) < 1.3 * (r_cov_i + r_cov_j).
   Ionic metals (Na,K,...) use near-zero cov_r to suppress spurious bonds.

2. TRANSITION METAL DATIVE BOND REMOVAL
   Keep only TM-C bonds (covalent sigma); remove all TM-X where X != C
   (dative bonds: the LP stays on the donor atom).

3. BOND ORDER UPGRADE  (NBO threshold-lowering, shortest-bond priority)
   Start all bonds as single. Greedily upgrade the SHORTEST eligible bond
   on the most deficient atom until all octets are satisfied.
   Deficiency uses plain (unformal) valence at this stage.

4. POST-UPGRADE FORMAL CHARGE REDISTRIBUTION  ← key NBO step
   After bond upgrade, hypervalent atoms (n_bonds > std_cap) carry a
   formal positive charge (excess = n_bonds - std_cap).  NBO maximises
   Lewis occupancy by moving this excess as formal -1 charges onto the
   most electronegative single-bonded neighbours.
   Two passes:
     A. Local:  for each hypervalent atom, distribute -1 to its directly
                bonded atoms that are (a) more electronegative, (b) have
                only a single bond to it (not a double/triple bond),
                (c) have bond_e < std_cap. Priority: O > N > S > C.
     B. Global: any unmatched positive charge (Na⁺, Pd⁺, I⁺ with no
                eligible local partner) is distributed to the globally
                most electronegative under-bonded atoms in the molecule.
   Each atom that receives a -1 formal charge gains one lone pair.
   TM lone pairs = (d_valence - n_C_bonds) // 2.

5. LONE PAIRS
   lp_i = max(0, eff_val_i - bond_e_i) // 2
   where eff_val = valence - max(0, n_bonds - std_cap) for hypervalent atoms,
         eff_val = valence + 1  for atoms receiving a formal -1 charge.

6. NON-LEWIS MODEL (empirical, fitted to nbo.log C8H10O)
   Valence non-Lewis = 0.017*n_sigma + 0.3263*n_pi + 0.020*n_lp
   Rydberg non-Lewis = 0.009*n_heavy + 0.003*n_H

Usage
-----
   python lewis_choose.py                  # demo: methylamine CH3NH2
   python lewis_choose.py mol.xyz          # XYZ file
   python lewis_choose.py mol.47           # NBO .47 file ($COORD block)
   python lewis_choose.py mol.xyz -1       # with molecular charge
"""

import sys, math
from collections import defaultdict

# ─── Periodic table ───────────────────────────────────────────────────────────

VALENCE = {
    'H':1,'He':2,
    'Li':1,'Be':2,'B':3,'C':4,'N':5,'O':6,'F':7,'Ne':8,
    'Na':1,'Mg':2,'Al':3,'Si':4,'P':5,'S':6,'Cl':7,'Ar':8,
    'K':1,'Ca':2,'Ga':3,'Ge':4,'As':5,'Se':6,'Br':7,'Kr':8,
    'Rb':1,'Sr':2,'In':3,'Sn':4,'Sb':5,'Te':6,'I':7,'Xe':8,
    'Cs':1,'Ba':2,'Tl':3,'Pb':4,'Bi':5,
    # TM: d+s electrons of neutral ground-state atom (used for LP count)
    'Sc':3,'Ti':4,'V':5,'Cr':6,'Mn':7,'Fe':8,'Co':9,'Ni':10,'Cu':11,'Zn':12,
    'Y':3,'Zr':4,'Nb':5,'Mo':6,'Tc':7,'Ru':8,'Rh':9,'Pd':10,'Ag':11,'Cd':12,
    'Hf':4,'Ta':5,'W':6,'Re':7,'Os':8,'Ir':9,'Pt':10,'Au':11,'Hg':12,
    'La':3,'Ce':4,'Pr':5,'Nd':6,'Pm':7,'Sm':8,'Eu':9,'Gd':10,'Tb':11,
    'Dy':12,'Ho':13,'Er':14,'Tm':15,'Yb':16,'Lu':17,
}

CORE_E = {
    'H':0,'He':0,
    'Li':2,'Be':2,'B':2,'C':2,'N':2,'O':2,'F':2,'Ne':2,
    'Na':10,'Mg':10,'Al':10,'Si':10,'P':10,'S':10,'Cl':10,'Ar':10,
    'K':18,'Ca':18,
    'Sc':18,'Ti':18,'V':18,'Cr':18,'Mn':18,'Fe':18,'Co':18,'Ni':18,'Cu':18,'Zn':18,
    'Ga':28,'Ge':28,'As':28,'Se':28,'Br':28,'Kr':28,
    'Rb':36,'Sr':36,
    'Y':36,'Zr':36,'Nb':36,'Mo':36,'Tc':36,'Ru':36,'Rh':36,'Pd':36,'Ag':36,'Cd':36,
    'In':46,'Sn':46,'Sb':46,'Te':46,'I':46,'Xe':46,
    'Cs':54,'Ba':54,
    'La':54,'Ce':54,'Pr':54,'Nd':54,'Pm':54,'Sm':54,'Eu':54,'Gd':54,'Tb':54,
    'Dy':54,'Ho':54,'Er':54,'Tm':54,'Yb':54,'Lu':54,
    'Hf':68,'Ta':68,'W':68,'Re':68,'Os':68,'Ir':68,'Pt':68,'Au':68,'Hg':68,
    'Tl':78,'Pb':78,'Bi':78,
}

COV_R = {
    'H':0.31,'He':0.28,
    'Li':0.25,'Be':0.96,'B':0.84,'C':0.76,'N':0.71,'O':0.66,'F':0.57,'Ne':0.58,
    'Na':0.25,'Mg':0.72,'Al':1.21,'Si':1.11,'P':1.07,'S':1.05,'Cl':1.02,'Ar':1.06,
    'K':0.25,'Ca':0.25,
    'Sc':1.44,'Ti':1.36,'V':1.25,'Cr':1.22,'Mn':1.19,'Fe':1.16,'Co':1.11,
    'Ni':1.10,'Cu':1.12,'Zn':1.18,
    'Ga':1.22,'Ge':1.20,'As':1.19,'Se':1.20,'Br':1.20,'Kr':1.16,
    'Rb':0.25,'Sr':0.25,
    'Y':1.62,'Zr':1.48,'Nb':1.37,'Mo':1.45,'Tc':1.56,'Ru':1.26,'Rh':1.35,
    'Pd':1.24,'Ag':1.45,'Cd':1.44,
    'In':1.42,'Sn':1.39,'Sb':1.39,'Te':1.38,'I':1.39,'Xe':1.40,
    'Cs':0.25,'Ba':0.25,
    'La':1.94,'Ce':1.83,'Pr':1.82,'Nd':1.81,'Pm':1.80,'Sm':1.80,'Eu':1.99,
    'Gd':1.79,'Tb':1.76,'Dy':1.75,'Ho':1.74,'Er':1.73,'Tm':1.72,'Yb':1.94,'Lu':1.72,
    'Hf':1.52,'Ta':1.46,'W':1.37,'Re':1.31,'Os':1.44,'Ir':1.41,
    'Pt':1.36,'Au':1.36,'Hg':1.32,
    'Tl':1.45,'Pb':1.46,'Bi':1.48,
}

# Standard covalent bond capacity (normal-valence, without hypervalency)
STD_CAP = {
    'H':1,'He':0,'Li':1,'Be':2,'B':3,'C':4,'N':3,'O':2,'F':1,'Ne':0,
    'Na':1,'Mg':2,'Al':3,'Si':4,'P':3,'S':2,'Cl':1,'Ar':0,
    'K':1,'Ca':2,'Ga':3,'Ge':4,'As':3,'Se':2,'Br':1,'Kr':0,
    'Rb':1,'Sr':2,'In':3,'Sn':4,'Sb':3,'Te':2,'I':1,'Xe':0,
    'Cs':1,'Ba':2,'Tl':3,'Pb':4,'Bi':3,
}

# Pauling electronegativity — determines priority for formal -1 charge assignment
ENEG = {
    'F':3.98,'O':3.44,'Cl':3.16,'N':3.04,'Br':2.96,'I':2.66,
    'S':2.58,'Se':2.55,'C':2.55,'P':2.19,'H':2.20,'As':2.18,
    'Te':2.10,'Si':1.90,'B':2.04,'Ge':2.01,'Sn':1.96,'Sb':2.05,
    'Pb':2.33,'Al':1.61,'Ga':1.81,'In':1.78,'Tl':2.04,
}

_TM = {
    'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
    'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
    'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
}

NL_SIGMA=0.017; NL_PI=0.3263; NL_LP=0.020; RYD_HEAVY=0.009; RYD_H=0.003


def is_TM(s):  return s in _TM
def is_H(s):   return s == 'H'
def octet(s):  return 0 if is_TM(s) else (2 if s in ('H','He','Li') else 8)
def std_cap(s):return STD_CAP.get(s, 4)
def eneg(s):   return ENEG.get(s, 2.0)


# ─── Geometry ─────────────────────────────────────────────────────────────────

def dist(a, b):
    return math.sqrt(sum((a[k]-b[k])**2 for k in range(3)))

def connectivity(atoms, coords, factor=1.3):
    bonds = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            ri = COV_R.get(atoms[i], 0.77)
            rj = COV_R.get(atoms[j], 0.77)
            # Use tighter factor for TM-C bonds to exclude long contacts
            # (e.g. amidinate C through N, or agostic interactions)
            if is_TM(atoms[i]) and atoms[j] == 'C' or is_TM(atoms[j]) and atoms[i] == 'C':
                f = 1.15
            else:
                f = factor
            if dist(coords[i], coords[j]) < f*(ri+rj):
                bonds.append((i, j))
    return bonds


# ─── Lewis structure engine ───────────────────────────────────────────────────

def _sigma_TM_C_set(atoms, bonds):
    """
    For each TM, identify which TM-C bonds are sigma (covalent) vs haptic/agostic (dative).

    A TM-C bond is treated as DATIVE (not sigma) when:
      (a) Haptic: C has another neighbour also bonded to the same TM
          (Cp, benzene, allyl coordination).
      (b) Agostic C-H: C has a H neighbour that is ALSO within bonding distance
          of the same TM AND C already has its full complement of non-TM bonds
          (bond_count_to_non_TM == STD_CAP[C]-1, i.e. C would be at 4 bonds
          if the TM bond were sigma). In this case the TM-C contact is agostic,
          driven by the C-H sigma donation, not a direct C-TM covalent bond.

    Only sigma bonds are kept (returned as a set of (tm, c) tuples).
    """
    # TM -> set of C atoms bonded to it
    tm_C = defaultdict(set)
    for i, j in bonds:
        if is_TM(atoms[i]) and atoms[j] == 'C':  tm_C[i].add(j)
        if is_TM(atoms[j]) and atoms[i] == 'C':  tm_C[j].add(i)
    # C adjacency (all bonds)
    c_adj = defaultdict(set)
    for i, j in bonds:
        c_adj[i].add(j); c_adj[j].add(i)
    # H adjacency — for agostic detection
    h_adj = defaultdict(set)   # h_idx -> set of atoms bonded to that H
    for i, j in bonds:
        if atoms[i] == 'H': h_adj[i].add(j)
        if atoms[j] == 'H': h_adj[j].add(i)
    # TM -> set of H atoms bonded to it (in raw connectivity)
    tm_H = defaultdict(set)
    for i, j in bonds:
        if is_TM(atoms[i]) and atoms[j] == 'H': tm_H[i].add(j)
        if is_TM(atoms[j]) and atoms[i] == 'H': tm_H[j].add(i)

    sigma = set()
    for tm, c_set in tm_C.items():
        for c in c_set:
            # Test (a): haptic — any C-neighbour of c is also TM-bonded
            haptic_nbrs = [k for k in c_adj[c]
                           if k != tm and not is_TM(atoms[k]) and k in c_set]
            if haptic_nbrs:
                continue   # haptic → dative, skip

            # Test (b): agostic C-H — c has a H neighbour that also bonds to tm
            # AND c's non-TM, non-H bonds already fill its standard capacity
            non_tm_nbrs = [k for k in c_adj[c] if k != tm and not is_TM(atoms[k])]
            h_nbrs_of_c = [k for k in non_tm_nbrs if atoms[k] == 'H']
            heavy_nbrs_of_c = [k for k in non_tm_nbrs if atoms[k] != 'H']

            # Check if any H on c is also bonded to tm (agostic H)
            agostic = False
            for h in h_nbrs_of_c:
                if h in tm_H[tm]:
                    agostic = True
                    break

            if agostic:
                # C is agostic donor via C-H → skip (treat as dative, not sigma)
                continue

            sigma.add((tm, c))
    return sigma


def remove_dative_tm_bonds(atoms, bonds):
    """
    Keep only covalent TM-C sigma bonds; drop:
      - TM-X bonds where X is not C (N, O, halide: dative)
      - haptic TM-C bonds (Cp, arene, allyl): dative from pi system
    """
    sigma = _sigma_TM_C_set(atoms, bonds)
    kept = []
    for i, j in bonds:
        si, sj = atoms[i], atoms[j]
        if is_TM(si) or is_TM(sj):
            tm, c = (i, j) if is_TM(si) else (j, i)
            if atoms[c] == 'C' and (tm, c) in sigma:
                kept.append((i, j))
            # else: dative (TM-X non-C, or haptic TM-C) → drop
        else:
            kept.append((i, j))
    return kept


def _deficiency_plain(i, atoms, bo, adj):
    """Deficiency using plain valence (before formal charge adjustment)."""
    if is_TM(atoms[i]): return 0
    val    = VALENCE.get(atoms[i], 4)
    bond_e = sum(bo[min(i,j), max(i,j)] for j in adj[i])
    lone_e = max(0, val - bond_e)
    needed = max(0, octet(atoms[i]) - 2*bond_e)
    return max(0, needed - lone_e)


def upgrade_bonds(atoms, coords, bo, adj):
    """Greedy bond-order upgrade: most deficient atom, shortest bond first."""
    improved = True
    while improved:
        improved = False
        cands = [i for i in range(len(atoms))
                 if not is_H(atoms[i]) and octet(atoms[i]) > 0]
        def _sort_key(i):
            d_val = -_deficiency_plain(i, atoms, bo, adj)
            # tie-break: atom whose SHORTEST eligible upgrade partner is shortest
            # (prioritises atoms connected via shorter bonds = higher pi density)
            nbr_dists = [
                dist(coords[i], coords[j])
                for j in adj[i]
                if not is_H(atoms[j]) and octet(atoms[j]) > 0
                and bo[min(i,j), max(i,j)] < 3
                and (
                    _deficiency_plain(j, atoms, bo, adj) > 0
                    or (sum(1 for k in adj[i] if not is_H(atoms[k])) == 1
                        and sum(1 for k in adj[j] if not is_H(atoms[k])) == 1
                        and sum(bo[min(j,k),max(j,k)] for k in adj[j]) < 3)
                )
            ]
            best_d = min(nbr_dists) if nbr_dists else math.inf
            return (d_val, best_d)
        heavy = sorted(cands, key=_sort_key)
        for i in heavy:
            if _deficiency_plain(i, atoms, bo, adj) == 0:
                break
            nbrs = sorted(
                [j for j in adj[i]
                 if not is_H(atoms[j])
                 and octet(atoms[j]) > 0
                 and bo[min(i,j), max(i,j)] < 3
                 and (
                     _deficiency_plain(j, atoms, bo, adj) > 0
                     # Both atoms must be terminal (1 non-H heavy neighbour each)
                     # to allow triple bonds in diatomics: CO, N2, CS, HCN, etc.
                     or (sum(1 for k in adj[i] if not is_H(atoms[k])) == 1
                         and sum(1 for k in adj[j] if not is_H(atoms[k])) == 1
                         and sum(bo[min(j,k),max(j,k)] for k in adj[j]) < 3)
                 )],
                key=lambda j: dist(coords[i], coords[j])
            )
            for j in nbrs:
                bo[min(i,j), max(i,j)] += 1
                improved = True
                break
            if improved:
                break


def redistribute_formal_charges(atoms, bo, adj, mol_charge):
    """
    Compute per-atom formal charges to maximise Lewis occupancy.

    Strategy
    --------
    Pass A  (local, hypervalent donors):
      Atoms with sum-of-bond-orders > std_cap are formally positive.
      Distribute -1 to their most electronegative single-bonded neighbours
      that have bond_e < std_cap (room for more electrons).

    Pass B  (global, molecular charge + residual):
      Any remaining imbalance (ionic metals, TM residuals, mol_charge) is
      distributed to the globally most-electronegative under-bonded atoms.
      Fallback: if std_cap-based eligibility finds nothing, use bond_e < VALENCE.

    Pass C  (LP transfer, non-hypervalent donors):
      When a less-electronegative atom still has electron density
      (eff_val - bond_e > 0) next to a more-electronegative atom with
      room (bond_e < std_cap), transfer one electron pair (fc donor +1,
      receiver -1).  Handles S in SO2/SO3, N-oxide, sulfonamide etc.

    Undo    (maximum-LP check):
      For each hypervalent atom, check whether keeping its positive fc
      increases total LP vs reverting to neutral. Undo if it doesn't help.
    """
    n = len(atoms)
    fc = [0] * n
    pos_remaining = {}

    # ── Identify hypervalent and ionic positive sources ───────────────────────
    for i in range(n):
        sym = atoms[i]
        if is_TM(sym) or is_H(sym): continue
        bond_e_i = sum(bo.get((min(i,k), max(i,k)), 1) for k in adj[i])
        cap = std_cap(sym)
        if bond_e_i > cap:
            exc = bond_e_i - cap
            fc[i] = +exc
            pos_remaining[i] = exc

    # Ionic metals: atoms with no Lewis bonds but non-zero valence.
    # EXCLUDE: atoms that appear to be X-type ligands (bonded to TM via dative bond
    # that was removed from Lewis). We detect these by checking the Lewis adj of TMs:
    # if a TM has adj entries, those are its sigma-bonded ligands. But we need the
    # REVERSE: atoms with no Lewis bonds that are TM neighbors in the molecule.
    # We use a simple heuristic: skip atoms with VALENCE ≤ 7 (halogens, main group)
    # that have no Lewis bonds — these should get fc=-1 from TM residual, not +VALENCE.
    # True ionic bare metals (Na⁺, K⁺, etc.) have very small cov_r and are handled
    # by the "near-zero cov_r" approach in connectivity.
    for i in range(n):
        sym = atoms[i]
        if not is_TM(sym) and not is_H(sym) and len(adj[i]) == 0:
            # Skip atoms that are likely TM ligands with dative bonds removed.
            # These are organics/halides: C, N, O, S, P, halogens.
            # True ionic metals (Li, Na, K, etc.) have VALENCE ≤ 2 and would
            # not normally be "ligands" in organometallic chemistry.
            v = VALENCE.get(sym, 0)
            if v > 0:
                # Only apply ionic positive charge for truly isolated main-group cations
                # (bare metal ions with no connectivity at all). Skip if it's a
                # non-metal atom that's likely a ligand (will get fc=-1 from TM residual).
                is_likely_ligand = sym in ('F', 'Cl', 'Br', 'I',
                                           'O', 'N', 'S', 'P', 'C',
                                           'Se', 'Te', 'As', 'Sb')
                if not is_likely_ligand:
                    fc[i] = +v
                    pos_remaining[i] = v

    # TM residual
    for i in range(n):
        sym = atoms[i]
        if is_TM(sym):
            val = VALENCE.get(sym, 0)
            nc  = sum(1 for j in adj[i] if atoms[j] == 'C')
            bond_e = nc
            lp_prelim = 0 if val <= 6 else max(0, val - nc) // 2
            residual = val - bond_e - 2*lp_prelim
            if residual > 0:
                fc[i] = +residual
                pos_remaining[i] = residual

    # ── Pass A: distribute hypervalent/ionic positives to local neighbours ────
    if pos_remaining:
        for src in sorted(pos_remaining.keys(),
                          key=lambda i: -VALENCE.get(atoms[i], 0)):
            if pos_remaining[src] <= 0: continue
            sym_s = atoms[src]
            en_s  = eneg(sym_s)
            nbr_cands = []
            for j in adj[src]:
                if fc[j] != 0: continue
                sym_j = atoms[j]
                if is_TM(sym_j) or is_H(sym_j): continue
                if eneg(sym_j) <= en_s: continue
                cap_j    = std_cap(sym_j)
                bond_e_j = sum(bo.get((min(j,k), max(j,k)), 1) for k in adj[j])
                if bond_e_j >= cap_j: continue
                key = (min(src,j), max(src,j))
                if bo.get(key, 1) > 1: continue
                nbr_cands.append((-eneg(sym_j), bond_e_j, j))
            nbr_cands.sort()
            for _, _, j in nbr_cands:
                if pos_remaining[src] <= 0: break
                fc[j] = -1
                pos_remaining[src] -= 1

    # ── Pass B: global redistribution for remaining imbalance ────────────────
    current_sum = sum(fc)
    needed_neg  = current_sum - mol_charge
    if needed_neg > 0:
        gcands = []
        for j in range(n):
            sym_j = atoms[j]
            if is_TM(sym_j) or is_H(sym_j): continue
            if fc[j] != 0: continue
            cap_j    = std_cap(sym_j)
            bond_e_j = sum(bo.get((min(j,k), max(j,k)), 1) for k in adj[j])
            if bond_e_j >= cap_j: continue
            gcands.append((-eneg(sym_j), bond_e_j, j))
        gcands.sort()
        assigned = 0
        for _, _, j in gcands:
            if assigned >= needed_neg: break
            fc[j] = -1; assigned += 1
        # Fallback: allow atoms at std_cap if still imbalanced
        if assigned < needed_neg:
            exp = []
            for j in range(n):
                sym_j = atoms[j]
                if is_TM(sym_j) or is_H(sym_j): continue
                if fc[j] != 0: continue
                val_j    = VALENCE.get(sym_j, 4)
                bond_e_j = sum(bo.get((min(j,k), max(j,k)), 1) for k in adj[j])
                if bond_e_j >= val_j: continue
                exp.append((-eneg(sym_j), -bond_e_j, j))
            exp.sort()
            for _, _, j in exp:
                if assigned >= needed_neg: break
                fc[j] = -1; assigned += 1

    # ── Pass C: LP transfer (non-hypervalent donors) ──────────────────────────
    # Handles S in SO2/SO3, N-oxides: atoms with EXCESS electron density (more
    # electrons than their normal bond pattern provides) next to more-electronegative
    # under-bonded neighbours.
    # IMPORTANT: Only trigger when the donor has electrons BEYOND its normal lone-pair
    # count, i.e. when bond_e_d < std_cap(donor). An atom at its normal bond count
    # (bond_e == std_cap) has exactly its expected lone pairs — nothing to transfer.
    # Example: N with 3 bonds (std_cap=3) has (5-3)=2 electrons = 1 normal LP.
    # This is NOT excess; Pass C must NOT fire for it (would wrongly give N fc=+1).
    changed = True
    while changed:
        changed = False
        for donor in range(n):
            sym_d = atoms[donor]
            if is_TM(sym_d) or is_H(sym_d): continue
            bond_e_d = sum(bo.get((min(donor,k), max(donor,k)), 1) for k in adj[donor])
            cap_d    = std_cap(sym_d)
            # Only fire when donor is UNDER-bonded (has excess electrons relative to
            # its normal valence pattern — not just its normal lone pairs)
            if bond_e_d >= cap_d: continue
            eff_v_d  = VALENCE.get(sym_d, 4) - fc[donor]
            remaining_e = max(0, eff_v_d - bond_e_d)
            if remaining_e == 0: continue
            en_d = eneg(sym_d)
            for recv in adj[donor]:
                sym_r = atoms[recv]
                if is_TM(sym_r) or is_H(sym_r): continue
                if fc[recv] != 0: continue
                if eneg(sym_r) <= en_d: continue
                cap_r    = std_cap(sym_r)
                bond_e_r = sum(bo.get((min(recv,k), max(recv,k)), 1) for k in adj[recv])
                if bond_e_r >= cap_r: continue
                key = (min(donor, recv), max(donor, recv))
                if bo.get(key, 1) > 1: continue
                fc[donor] += 1
                fc[recv]   = -1
                changed = True
                break
            if changed: break

    # ── Undo / Reassign: maximise total LP for hypervalent atoms ─────────────
    # For each atom with fc > 0, test three alternatives:
    #   (a) keep current fc (baseline)
    #   (b) zero out fc and matched receivers (revert to neutral)
    #   (c) flip fc to -1 and zero matched receivers (makes atom a net acceptor)
    # Choose the alternative with highest total LP that preserves charge balance.
    def total_lp_with_fc(fc_trial):
        total = 0
        for i, sym in enumerate(atoms):
            if is_TM(sym): continue
            ev = VALENCE.get(sym, 4) - fc_trial[i]
            be = sum(bo.get((min(i,k), max(i,k)), 1) for k in adj[i])
            total += max(0, ev - be) // 2
        return total

    for i in range(n):
        if fc[i] <= 0 or is_TM(atoms[i]): continue
        excess = fc[i]
        current_lp = total_lp_with_fc(fc)

        # Build trial (b): zero out atom i and its matched -1 receivers
        trial_b = fc[:]
        trial_b[i] = 0
        undone = 0
        for j in list(adj[i]) + list(range(n)):
            if undone >= excess: break
            if trial_b[j] == -1 and j != i:
                trial_b[j] = 0; undone += 1

        # Build trial (c): flip atom i to -1 and zero all matched receivers
        trial_c = fc[:]
        trial_c[i] = -1
        undone = 0
        for j in list(adj[i]) + list(range(n)):
            if undone >= excess + 1: break   # need to undo excess+1 negatives
            if trial_c[j] == -1 and j != i:
                trial_c[j] = 0; undone += 1

        best = fc
        for trial in (trial_b, trial_c):
            if sum(trial) != mol_charge: continue
            if total_lp_with_fc(trial) > total_lp_with_fc(best):
                best = trial
        if best is not fc:
            fc = best

    return fc


def compute_lone_pairs(atoms, bo, adj, fc):
    lp = {}
    for i, sym in enumerate(atoms):
        if is_TM(sym):
            val = VALENCE.get(sym, 0)
            # After dative bond removal, adj[i] contains only sigma TM-C bonds
            nc  = sum(1 for j in adj[i] if atoms[j] == 'C')
            # Early TMs (val <= 6) donate all valence electrons → 0 LP
            lp[i] = 0 if val <= 6 else max(0, val - nc) // 2
        else:
            val    = VALENCE.get(sym, 4)
            eff_v  = val - fc[i]   # fc[i]=+exc reduces eff_val; fc[i]=-1 increases it
            bond_e = sum(bo[min(i,j), max(i,j)] for j in adj[i])
            lp[i]  = max(0, eff_v - bond_e) // 2
    return lp


def _is_stuck(i, atoms, bo, adj):
    """True if atom i is deficient but has no eligible upgrade partner."""
    if is_H(atoms[i]) or octet(atoms[i]) == 0: return False
    if _deficiency_plain(i, atoms, bo, adj) == 0: return False
    for j in adj[i]:
        if is_H(atoms[j]) or octet(atoms[j]) == 0: continue
        if bo[min(i,j), max(i,j)] >= 3: continue
        if (_deficiency_plain(j, atoms, bo, adj) > 0
                or (sum(1 for k in adj[i] if not is_H(atoms[k])) == 1
                    and sum(1 for k in adj[j] if not is_H(atoms[k])) == 1
                    and sum(bo[min(j,k),max(j,k)] for k in adj[j]) < 3)):
            return False   # has an eligible partner
    return True  # no eligible partner → stuck


def _ring_atoms(adj, atoms):
    """
    Return set of atom indices that participate in at least one ring (cycle).
    Uses iterative DFS to avoid Python recursion limits.
    """
    n = len(atoms)
    in_ring = set()
    visited = {}   # node -> int (discovery order)
    parent  = {}   # node -> parent node
    on_stack = set()
    counter = [0]

    for start in range(n):
        if is_H(atoms[start]) or start in visited:
            continue
        # Iterative DFS from start
        # stack entries: (node, iterator-over-neighbours)
        stack = [(start, iter(adj[start]))]
        visited[start] = counter[0]; counter[0] += 1
        on_stack.add(start)

        while stack:
            v, children = stack[-1]
            advanced = False
            for w in children:
                if is_H(atoms[w]):
                    continue
                if w not in visited:
                    parent[w] = v
                    visited[w] = counter[0]; counter[0] += 1
                    on_stack.add(w)
                    stack.append((w, iter(adj[w])))
                    advanced = True
                    break
                elif w != parent.get(v, -1) and w in on_stack:
                    # Back edge v→w: trace the cycle path and mark all atoms
                    in_ring.add(v); in_ring.add(w)
                    x = v
                    while x != w and x != -1:
                        in_ring.add(x)
                        x = parent.get(x, -1)
            if not advanced:
                on_stack.discard(v)
                stack.pop()

    return in_ring


def _kekulé_swap(atoms, coords, bo, adj):
    """
    Fix stuck ring atoms after the greedy upgrade by rotating Kekulé doubles.

    When two non-alternating double bonds are placed in a ring (e.g. C1=C2 and
    C4=C5 leaving C3 and C6 with no upgrade partner), swap one double bond to
    open a path, then re-run the greedy upgrade.

    Only atoms that are part of a ring (detected by DFS cycle-finding) are
    eligible for swapping.  Terminal atoms (carboxylate O, amidinate N, etc.)
    and acyclic chain atoms are left for formal-charge redistribution.
    """
    ring_set = _ring_atoms(adj, atoms)
    MAX_SWAPS = len(atoms)   # safety cap
    for _ in range(MAX_SWAPS):
        # Find stuck ring atoms only
        stuck = [i for i in range(len(atoms))
                 if i in ring_set
                 and _is_stuck(i, atoms, bo, adj)
                 and not is_H(atoms[i])
                 and octet(atoms[i]) > 0]
        if not stuck:
            break

        n_stuck = len(stuck)
        swapped = False
        for i in stuck:
            if swapped: break
            for j in adj[i]:
                if is_H(atoms[j]) or swapped: continue
                key_ij = (min(i, j), max(i, j))
                if bo[key_ij] >= 2: continue   # already doubled
                # Look for a double bond on j to a third atom k
                for k in adj[j]:
                    if k == i: continue
                    key_jk = (min(j, k), max(j, k))
                    if bo[key_jk] != 2: continue
                    # Trial swap: remove D j-k, add D i-j
                    bo_trial = dict(bo)
                    bo_trial[key_jk] = 1
                    bo_trial[key_ij] = 2
                    # Count stuck ring atoms after swap + re-upgrade
                    bo_copy = dict(bo_trial)
                    upgrade_bonds(atoms, coords, bo_copy, adj)
                    n_trial = sum(1 for x in ring_set
                                  if _is_stuck(x, atoms, bo_copy, adj)
                                  and not is_H(atoms[x])
                                  and octet(atoms[x]) > 0)
                    if n_trial < n_stuck:
                        # Accept: apply swap and re-upgrade the real bo dict
                        bo[key_jk] = 1
                        bo[key_ij] = 2
                        upgrade_bonds(atoms, coords, bo, adj)
                        swapped = True
                        break
                if swapped: break
        if not swapped:
            break   # no beneficial swap found, stop


def _shift_acyl_lp(atoms, bo, adj, lp, fc):
    """
    Post-processing: move lone pairs from carbon to oxygen in acyl/enolate systems.

    Problem: the greedy bond-upgrade algorithm sometimes builds both C=O double bonds
    in a symmetric acyl fragment (beta-diketonate, malonyl, carboxylate) and assigns
    the resulting lone pair to the central carbon instead of the more electronegative
    oxygen.  For example:
        O=C-CH-C=O   (LP on CH)   ← wrong
        O=C-CH=C-O⁻  (LP on O⁻)  ← correct

    Detection:
        Carbon i has LP ≥ 1  AND  has no own double bond in bo
        AND at least one neighbour j is a carbonyl carbon (C bonded to O with bo ≥ 2)
        AND that oxygen has bo ≥ 2 to j (is the double-bond O)

    Action (first eligible pair found; iterates until no more):
        1. Decrease bo[j, O] from 2 → 1  (break C=O double bond)
        2. Increase bo[i, j] from 1 → 2  (make C=C double bond)
        3. Transfer LP: lp[i] -= 1, lp[O] += 1
        4. Update fc: O gets fc = -1 (now anionic with 1 bond), i gets fc restored.

    Only applies to carbon; does not touch TMs or heteroatoms.
    """
    changed = True
    while changed:
        changed = False
        for i, sym_i in enumerate(atoms):
            if sym_i != 'C':
                continue
            if lp.get(i, 0) < 1:
                continue
            # i must have no existing double bond
            if any(bo.get((min(i, k), max(i, k)), 0) >= 2 for k in adj[i]):
                continue
            # Look for a neighbour j that is a carbonyl C (C bonded to O via double bond)
            for j in adj[i]:
                if atoms[j] != 'C':
                    continue
                key_ij = (min(i, j), max(i, j))
                if bo.get(key_ij, 0) != 1:  # must be a single bond currently
                    continue
                # Find the carbonyl oxygen on j
                o_cand = None
                for k in adj[j]:
                    if atoms[k] == 'O':
                        key_jk = (min(j, k), max(j, k))
                        if bo.get(key_jk, 0) >= 2:
                            o_cand = k
                            break
                if o_cand is None:
                    continue
                # Check that the oxygen o_cand has no other bonds (it is terminal)
                other_bonds = [k for k in adj[o_cand] if k != j]
                if other_bonds:
                    continue  # O has other substituents, don't touch it
                # Perform the shift:
                key_jO = (min(j, o_cand), max(j, o_cand))
                bo[key_jO]  -= 1   # C=O → C-O
                bo[key_ij]  += 1   # C-C → C=C
                lp[i]       -= 1   # remove LP from carbon
                lp[o_cand]  += 1   # add LP to oxygen
                fc[o_cand]   = -1  # O now anionic (1 bond, gained LP)
                # Recalculate carbon i's fc (now it has one more bond)
                fc[i] = max(0, fc[i])  # i is now doubly bonded, generally neutral
                changed = True
                break
            if changed:
                break


def _find_simple_rings(adj, atoms, max_size=7):
    """
    Find all simple rings up to max_size using iterative DFS.
    Returns list of frozensets of atom indices (one per unique ring).
    Only considers non-H atoms.
    """
    n = len(atoms)
    rings = []
    seen_rings = set()

    for start in range(n):
        if is_H(atoms[start]):
            continue
        # DFS with path tracking
        stack = [(start, [start], {start})]
        while stack:
            node, path, visited = stack.pop()
            for nb in adj[node]:
                if is_H(atoms[nb]):
                    continue
                if nb == start and len(path) >= 4:
                    key = frozenset(path)
                    if key not in seen_rings:
                        seen_rings.add(key)
                        rings.append(list(path))
                elif nb not in visited and len(path) < max_size:
                    stack.append((nb, path + [nb], visited | {nb}))
    return rings


def _is_sp2_candidate(i, atoms, adj):
    """
    Return True if atom i could be sp2 (part of aromatic/conjugated system).
    Criteria: C, N, O, S with ≤ 3 heavy neighbors and octet-requiring.
    """
    sym = atoms[i]
    if sym not in ('C', 'N', 'O', 'S'):
        return False
    heavy_nbrs = [j for j in adj[i] if not is_H(atoms[j])]
    # C with 3 heavy neighbors: sp2 candidate
    # N with 2-3 heavy neighbors: sp2 candidate
    # O with 2 heavy neighbors that are aromatic: sp2 (furan O)
    if sym == 'C' and len(heavy_nbrs) <= 3:
        return True
    if sym == 'N' and 2 <= len(heavy_nbrs) <= 3:
        return True
    if sym == 'O' and len(heavy_nbrs) == 2:
        return True
    if sym == 'S' and len(heavy_nbrs) == 2:
        return True
    return False


def _enforce_aromatic_rings(atoms, coords, bo, adj):
    """
    Detect and fix aromatic/conjugated rings with wrong Kekulé patterns.

    Strategy: Find all 5- and 6-membered sp2 rings, group fused systems together,
    and jointly optimise all ring bonds within each fused system by exhaustive
    enumeration of valid alternating-double patterns across all rings simultaneously.
    For small systems (≤ 12 ring bonds) this is fast; larger systems use greedy iteration.
    """
    def bo_ij(i, j):
        return bo.get((min(i, j), max(i, j)), 0)

    def total_be(i):
        return sum(bo_ij(i, k) for k in adj[i])

    def target_be(sym):
        return 4 if sym == 'C' else (3 if sym in ('N', 'P') else 2)

    def atom_ok(i):
        return total_be(i) == target_be(atoms[i])

    def atom_deficient(i):
        sym = atoms[i]
        if is_H(sym) or is_TM(sym):
            return False
        return total_be(i) < target_be(sym)

    # Find candidate aromatic rings
    all_rings = _find_simple_rings(adj, atoms, max_size=7)
    aromatic_rings = []
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        if not all(_is_sp2_candidate(i, atoms, adj) for i in ring):
            continue
        cn_count = sum(1 for i in ring if atoms[i] in ('C', 'N'))
        if cn_count < 4:
            continue
        aromatic_rings.append(frozenset(ring))

    # Deduplicate
    aromatic_rings = list(dict.fromkeys(aromatic_rings))

    if not aromatic_rings:
        return

    # Group into fused systems: two rings are fused if they share ≥ 2 atoms (an edge)
    def find_fused_groups(rings):
        n = len(rings)
        group = list(range(n))
        def find(x):
            while group[x] != x:
                group[x] = group[group[x]]
                x = group[x]
            return x
        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                group[rx] = ry
        for i in range(n):
            for j in range(i + 1, n):
                if len(rings[i] & rings[j]) >= 2:
                    union(i, j)
        groups = defaultdict(list)
        for i in range(n):
            groups[find(i)].append(i)
        return list(groups.values())

    fused_groups = find_fused_groups(aromatic_rings)

    for group_indices in fused_groups:
        group_rings = [aromatic_rings[i] for i in group_indices]

        # Collect all atoms and all ring bonds in this fused system
        all_atoms_in_group = set()
        for r in group_rings:
            all_atoms_in_group |= r

        # Check if any atom in this group is deficient
        if not any(atom_deficient(i) for i in all_atoms_in_group):
            continue

        # Collect all bonds that are part of ANY ring in this group
        # We need to enumerate double-bond placement on these bonds.
        # Bond: internal (both ends in group) vs exo (one end outside group).
        group_bonds = set()
        for r in group_rings:
            r_list = sorted(r)
            # Find ring edges by connectivity
            for i in r_list:
                for j in adj[i]:
                    if j in r and (min(i, j), max(i, j)) not in group_bonds:
                        # Check if this edge is inside any ring
                        for ring in group_rings:
                            if i in ring and j in ring:
                                group_bonds.add((min(i, j), max(i, j)))
                                break

        group_bonds = sorted(group_bonds)
        n_bonds = len(group_bonds)

        # Exo bond-electrons for each atom in group
        exo_be = {}
        for i in all_atoms_in_group:
            exo_be[i] = sum(bo_ij(i, k) for k in adj[i] if k not in all_atoms_in_group)

        # Identify exocyclic double bonds that could be relocated to allow ring aromaticity
        # Two types:
        #   1. Terminal C=O/N bonds (original logic)
        #   2. C=C bonds where the exo C has alternative double-bond partners (NHC junction)
        exo_double_o = {}
        for i in all_atoms_in_group:
            if atoms[i] != 'C':
                continue
            for k in adj[i]:
                if k in all_atoms_in_group:
                    continue
                key = (min(i, k), max(i, k))
                if bo.get(key, 1) >= 2:
                    other = [m for m in adj[k] if m != i and not is_H(atoms[m])]
                    other_all = [m for m in adj[k] if m != i]
                    if atoms[k] in ('O', 'N', 'S') and len(other) == 0:
                        # Terminal O=C or N=C (no heavy substituents) → candidate
                        exo_double_o[i] = (k, bo[key])
                        break
                    elif atoms[k] in ('N',) and len(other) <= 1:
                        # N with at most one heavy substituent (e.g. N7-N6 where N6
                        # has no further substituents). The exo C=N can be relocated
                        # to enable ring aromaticity; the N7-N6 fragment just shifts.
                        exo_double_o[i] = (k, bo[key])
                        break
                    elif atoms[k] == 'C' and len(other) >= 1:
                        # Exo C=C where exo-C has other sp2 neighbors it can use
                        exo_nbrs_nonring = [m for m in other if m not in all_atoms_in_group
                                            and atoms[m] == 'C']
                        if exo_nbrs_nonring:
                            exo_double_o[i] = (k, bo[key])
                            break

        def score_full(bond_orders_dict, exo_brk):
            """Score a complete assignment of bond orders to group_bonds."""
            score = 0
            for i in all_atoms_in_group:
                sym = atoms[i]
                tgt = target_be(sym)
                grp_be = sum(bond_orders_dict.get((min(i,k),max(i,k)), 1)
                             for k in adj[i] if k in all_atoms_in_group)
                my_exo = exo_be[i]
                if i in exo_brk and i in exo_double_o:
                    my_exo -= 1
                total = grp_be + my_exo
                was_def = atom_deficient(i)
                if total == tgt:
                    score += 3 if was_def else 1
                elif total > tgt:
                    score -= 2
            for ri in exo_brk:
                if ri in exo_double_o:
                    score += 1
            return score

        # Current state
        cur_bo = {b: bo.get(b, 1) for b in group_bonds}
        cur_mask_val = score_full(cur_bo, frozenset())
        best_score = cur_mask_val
        best_bo = dict(cur_bo)
        best_exo = frozenset()

        # For small systems enumerate exhaustively; for larger use greedy
        candidate_exo_sets = [frozenset()] + [frozenset([ri]) for ri in exo_double_o]

        if n_bonds <= 14:
            # Exhaustive: try all 2^n_bonds patterns
            for exo_brk in candidate_exo_sets:
                for mask in range(1 << n_bonds):
                    trial = {}
                    for idx, b in enumerate(group_bonds):
                        trial[b] = 2 if (mask >> idx) & 1 else 1
                    # Validity: no atom in group exceeds target
                    valid = True
                    for i in all_atoms_in_group:
                        sym = atoms[i]
                        tgt = target_be(sym)
                        grp_be = sum(trial.get((min(i,k),max(i,k)),1)
                                     for k in adj[i] if k in all_atoms_in_group)
                        my_exo = exo_be[i]
                        if i in exo_brk and i in exo_double_o:
                            my_exo -= 1
                        if grp_be + my_exo > tgt + 1:  # allow slight over for scoring to penalise
                            valid = False; break
                    if not valid:
                        continue
                    s = score_full(trial, exo_brk)
                    if s > best_score:
                        best_score = s
                        best_bo = dict(trial)
                        best_exo = exo_brk
        else:
            # Greedy: iterate ring-by-ring within the group
            for _iter in range(n_bonds):
                changed = False
                for ring in group_rings:
                    if not any(atom_deficient(i) for i in ring):
                        continue
                    r_list = sorted(ring)
                    r_edges = []
                    for i in r_list:
                        for j in adj[i]:
                            if j in ring and j > i:
                                r_edges.append((i, j))
                    n_re = len(r_edges)
                    best_local = score_full({b: bo.get(b,1) for b in group_bonds}, frozenset())
                    best_local_bo = None
                    for mask in range(1 << n_re):
                        trial = dict(cur_bo)
                        for idx, (a, b) in enumerate(r_edges):
                            trial[(min(a,b),max(a,b))] = 2 if (mask>>idx)&1 else 1
                        s = score_full(trial, frozenset())
                        if s > best_local:
                            best_local = s; best_local_bo = dict(trial)
                    if best_local_bo:
                        cur_bo = best_local_bo
                        if best_local > best_score:
                            best_score = best_local
                            best_bo = dict(cur_bo)
                        changed = True
                if not changed:
                    break

        # Apply best if improvement
        if best_score > cur_mask_val:
            for key, order in best_bo.items():
                bo[key] = order
            for ri in best_exo:
                if ri in exo_double_o:
                    exo_idx, _ = exo_double_o[ri]
                    key = (min(ri, exo_idx), max(ri, exo_idx))
                    if bo.get(key, 1) >= 2:
                        bo[key] -= 1


def check_octet_violations(atoms, bo, adj, lp=None, fc=None):
    """
    Check the Lewis structure for organic atoms that violate the octet rule.
    An atom violates if bond_e + 2*LP ≠ VALENCE[sym] - fc[i] (expected valence electrons).

    Exceptions:
      - TM atoms (variable electron count)
      - H atoms (duet rule)
      - Atoms where the neutral electron count (bond_e + 2*lp) already matches VALENCE,
        but fc != 0 was assigned as part of ionic TM charge redistribution — these are
        internally consistent Lewis structures where the formal charge is a bookkeeping
        artefact of distributing TM's oxidation state, not a true violation.
    """
    violations = []
    check_syms = {'C', 'N', 'O', 'S', 'P', 'B', 'Si', 'Se', 'Te', 'F', 'Cl', 'Br', 'I'}

    for i, sym in enumerate(atoms):
        if sym not in check_syms:
            continue
        bond_e = sum(bo.get((min(i, k), max(i, k)), 0) for k in adj[i])
        lp_cnt = lp.get(i, 0) if lp else 0
        fc_val = fc[i] if fc else 0
        total_e = bond_e + 2 * lp_cnt
        # Expected valence electrons = VALENCE adjusted for formal charge
        exp_val = VALENCE.get(sym, 4) - fc_val

        if total_e == exp_val:
            continue  # no violation

        # Special case: if the atom has fc != 0 but the NEUTRAL count is already
        # consistent (total_e == VALENCE), the fc is a bookkeeping artefact from
        # TM charge redistribution — not a true valence violation.
        neutral_val = VALENCE.get(sym, 4)
        if fc_val != 0 and total_e == neutral_val:
            continue  # neutral count OK; fc is TM ionic bookkeeping, not a violation

        # For atoms that can be trivalent/hypervalent (Si, P, Ge, Sn, B, etc.):
        # only flag if bond_e < 3 (truly undercoordinated, not just hypervalent).
        if sym in ('B', 'Si', 'P', 'Ge', 'Sn', 'As', 'Sb', 'Se', 'Te', 'Pb', 'Bi'):
            if bond_e >= 3:
                continue   # 3+ bonds is acceptable for these elements as ligands

        violations.append((i, f"{sym}{i+1}", bond_e, lp_cnt, exp_val, fc_val))
    return violations


def find_best_lewis(atoms, coords, charge=0):
    raw_bonds = connectivity(atoms, coords)
    if not raw_bonds:
        # Single atom or isolated ion: no bonds, compute LP from total electron count
        core_e  = sum(CORE_E.get(a, 0) for a in atoms)
        val_e   = sum(VALENCE.get(a, 0) for a in atoms) - charge
        total_e = core_e + val_e
        lp = {}
        for i, sym in enumerate(atoms):
            val = VALENCE.get(sym, 0) - charge  # all charge goes to this atom
            lp[i] = max(0, val) // 2
        fc = [0] * len(atoms)
        stats = dict(total_e=float(total_e), lewis_e=float(total_e),
                     val_nl=0.0, ryd_nl=0.0, charge=float(charge))
        return {}, lp, fc, stats

    bonds = remove_dative_tm_bonds(atoms, raw_bonds)

    adj = defaultdict(list)
    for i, j in bonds:
        adj[i].append(j); adj[j].append(i)

    bo = {(min(i,j), max(i,j)): 1 for i,j in bonds}

    # Step 1: bond upgrade (plain valence, shortest-bond rule)
    upgrade_bonds(atoms, coords, bo, adj)

    # Step 1b: Kekulé rotation — fix stuck ring atoms
    # When a ring carbon is still deficient after the greedy upgrade
    # (no eligible upgrade partners), the greedy placed non-alternating doubles.
    # Swap one conflicting double to free up a path, then re-upgrade.
    # Only applies to non-terminal heavy atoms (>=2 non-H heavy neighbours),
    # i.e. ring atoms. Terminal atoms (=1 non-H heavy nbr) like carboxylate O
    # are correctly left for the formal-charge redistribution in Step 2.
    _kekulé_swap(atoms, coords, bo, adj)

    # Step 1c: Aromatic ring enforcement
    # The Kekulé swap + greedy can still leave aromatic rings with wrong Kekulé
    # patterns when ring fusions or exocyclic double bonds conflict.
    # Enumerate all alternating-double patterns per candidate ring; pick the best.
    _enforce_aromatic_rings(atoms, coords, bo, adj)

    # Step 1d: CO/CS triple bond enforcement for metal carbonyls.
    # The greedy upgrade leaves M–C=O (double C=O) instead of M–C≡O (triple).
    # When a C bonded to a TM has a terminal O/S at bo=2 and is still deficient
    # (bond_e < 4), upgrade C=O → C≡O to give the correct carbonyl Lewis structure.
    raw_bonds_for_co = connectivity(atoms, coords)
    adj_raw_co = defaultdict(set)
    for a, b in raw_bonds_for_co:
        adj_raw_co[a].add(b); adj_raw_co[b].add(a)
    for c_idx, sym_c in enumerate(atoms):
        if sym_c != 'C':
            continue
        bond_e_c = sum(bo.get((min(c_idx, k), max(c_idx, k)), 0) for k in adj[c_idx])
        if bond_e_c >= 4:
            continue
        tm_nbrs = [k for k in adj_raw_co[c_idx] if is_TM(atoms[k])]
        if not tm_nbrs:
            continue
        for o_idx in adj[c_idx]:
            sym_o = atoms[o_idx]
            if sym_o not in ('O', 'S', 'N'):
                continue
            key_co = (min(c_idx, o_idx), max(c_idx, o_idx))
            if bo.get(key_co, 0) != 2:
                continue
            o_heavy = [k for k in adj[o_idx] if k != c_idx and not is_H(atoms[k])]
            if o_heavy:
                continue  # not terminal
            bo[key_co] = 3  # C=O → C≡O
            break

    # Step 1e-pre: Fix nitrate/nitrite N=O double bond.
    # When N is bonded exclusively to O (all neighbours are O) and has no double bond,
    # force one N=O upgrade so the Lewis structure has the correct N=O double bond.
    # This handles NO3-, NO2-, and similar oxyanion ligands where the greedy algorithm
    # distributed all lone pairs onto O atoms, leaving N with only single bonds.
    for n_idx, sym in enumerate(atoms):
        if sym != 'N':
            continue
        n_nbrs = list(adj[n_idx])
        if not n_nbrs:
            continue
        # N bonded exclusively to O (no C, no H, no other heteroatoms)
        if not all(atoms[k] == 'O' for k in n_nbrs):
            continue
        bond_e_n = sum(bo.get((min(n_idx, k), max(n_idx, k)), 0) for k in n_nbrs)
        has_double = any(bo.get((min(n_idx, k), max(n_idx, k)), 0) >= 2 for k in n_nbrs)
        if has_double:
            continue   # already has N=O, fine
        # N has only single bonds to O — upgrade the shortest N-O to a double bond
        best_o = min(n_nbrs, key=lambda k: dist(coords[n_idx], coords[k]))
        key_no = (min(n_idx, best_o), max(n_idx, best_o))
        bo[key_no] = 2

    # Step 1e-azide: Fix azide (N3-) N–N–N bond orders.
    # The greedy upgrade leaves linear N-N-N chains with a 2+1 pattern (double+single)
    # instead of the correct 2+2 (both double bonds). This happens because the central
    # N already has bond_e=3 (one double + one single), which satisfies its lone-pair
    # deficiency criterion (lone_e=2 fills the gap), so the greedy considers it done.
    #
    # Correct azide Lewis (most common resonance):
    #   Nα=Nβ+=Nγ-  :  bo(α-β)=2, bo(β-γ)=2, fc(β)=+1, fc(γ)=-1
    #   This gives Nβ bond_e=4, exactly filling its octet with no lone pairs (fc=+1).
    #
    # Detection: terminal N (1 heavy neighbour) bonded to a central N (2 N neighbours)
    #   bonded to another terminal N — all N, no other heavy atoms. Linear N-N-N chain.
    # Action: force both N-N bonds to bo=2.
    for n_cen in range(len(atoms)):
        if atoms[n_cen] != 'N':
            continue
        n_nbrs_cen = [k for k in adj[n_cen] if not is_H(atoms[k])]
        if len(n_nbrs_cen) != 2:
            continue
        # Central N must have exactly two heavy neighbours, both N
        if not all(atoms[k] == 'N' for k in n_nbrs_cen):
            continue
        n_a, n_b = n_nbrs_cen
        # Each terminal N must have only one heavy neighbour (the central N)
        n_a_heavy = [k for k in adj[n_a] if not is_H(atoms[k])]
        n_b_heavy = [k for k in adj[n_b] if not is_H(atoms[k])]
        if len(n_a_heavy) != 1 or len(n_b_heavy) != 1:
            continue
        # This is a linear N-N-N chain (azide-type) — force both bonds to bo=2
        key_ca = (min(n_cen, n_a), max(n_cen, n_a))
        key_cb = (min(n_cen, n_b), max(n_cen, n_b))
        bo[key_ca] = 2
        bo[key_cb] = 2

    # Step 1e: N–C–N amidinate / guanidinate double-bond shift (single pass).
    # In N-C-N chelate fragments the greedy algorithm places double bonds on BOTH
    # terminal N atoms (N=C…N…C=N), leaving the central N with only 2 bonds.
    # Correct: one terminal N keeps the double bond; the other loses it to the
    # central N, making the central N fully bonded and the terminal N X-type.
    # GUARD: only apply when the central N is NOT part of any ring.
    # Run a SINGLE PASS (no while loop) to avoid oscillation / infinite loops.
    ring_set_amidinate = _ring_atoms(adj, atoms)
    shifted_bonds = set()   # track bonds already shifted to prevent reversal

    for n_cen in range(len(atoms)):
        if atoms[n_cen] != 'N':
            continue
        if n_cen in ring_set_amidinate:
            continue
        bond_e_n = sum(bo.get((min(n_cen, k), max(n_cen, k)), 0) for k in adj[n_cen])
        if bond_e_n != 2:
            continue
        c_nbrs = [k for k in adj[n_cen] if atoms[k] == 'C']
        if len(c_nbrs) < 1:
            continue
        # Find a C neighbor that has a double bond to a terminal N (not n_cen)
        shifted = False
        for c1 in c_nbrs:
            if shifted:
                break
            for n_term in adj[c1]:
                if n_term == n_cen or atoms[n_term] != 'N':
                    continue
                if n_term in ring_set_amidinate:
                    continue
                key_c1_nt = (min(c1, n_term), max(c1, n_term))
                key_c1_nc = (min(c1, n_cen),  max(c1, n_cen))
                if (key_c1_nt in shifted_bonds or key_c1_nc in shifted_bonds):
                    continue   # already shifted this bond — skip to prevent oscillation
                if bo.get(key_c1_nt, 0) >= 2 and bo.get(key_c1_nc, 0) == 1:
                    bo[key_c1_nt] -= 1
                    bo[key_c1_nc] += 1
                    shifted_bonds.add(key_c1_nt)
                    shifted_bonds.add(key_c1_nc)
                    shifted = True
                    break

    # Step 1f: Dihydrogen (η²-H₂) H-H bond addition.
    # Two H atoms both coordinated to the same TM with no other Lewis bonds, and within
    # ~1.15 Å of each other, form a dihydrogen ligand. Their H-H σ-bond must appear in
    # the $CHOOSE BOND block even though the normal connectivity threshold (1.3×r_cov)
    # misses the stretched H-H bond (typically 0.85–1.05 Å in complexes vs 0.74 Å in free H₂).
    # Build TM-H adjacency from the raw (pre-dative-removal) connectivity for detection.
    _raw_bonds_hh = connectivity(atoms, coords)
    _adj_raw_hh = defaultdict(set)
    for _a, _b in _raw_bonds_hh:
        _adj_raw_hh[_a].add(_b); _adj_raw_hh[_b].add(_a)
    for _tm in range(len(atoms)):
        if not is_TM(atoms[_tm]):
            continue
        _h_candidates = [_h for _h in _adj_raw_hh[_tm] if atoms[_h] == 'H'
                         and len(adj[_h]) == 0]   # H with no Lewis bonds
        for _i, _hi in enumerate(_h_candidates):
            for _hj in _h_candidates[_i+1:]:
                _d = dist(coords[_hi], coords[_hj])
                if _d < 1.15:   # dihydrogen H-H bond threshold
                    _key = (min(_hi, _hj), max(_hi, _hj))
                    if _key not in bo:
                        bo[_key] = 1
                        adj[_hi].append(_hj); adj[_hj].append(_hi)

    # Step 1g: Azide (N3) bond-order correction.
    #
    # Free azide N3^- is isoelectronic with CO2: the correct Lewis structure is
    # [:N=N+=N:]^- (two cumulated double bonds, middle N positive, both ends negative).
    # When coordinated end-on to a metal through N_alpha, the dative Zn-N bond is
    # removed from the Lewis structure, leaving N_alpha with only one Lewis bond
    # (to N_beta). The greedy upgrade then incorrectly places the double bond on the
    # DISTAL N_beta=N_gamma pair (because N_gamma appears maximally deficient first),
    # giving the wrong pattern N_alpha-N_beta=N_gamma instead of the correct
    # N_alpha=N_beta=N_gamma.
    #
    # Detection: a linear N-N-N chain where
    #   (a) N_alpha has exactly ONE Lewis bond (to N_beta only) — TM bond was removed
    #   (b) N_beta has exactly TWO Lewis bonds (N_alpha and N_gamma)
    #   (c) N_gamma has exactly ONE Lewis bond (to N_beta only) — terminal
    #   (d) N_alpha is within bonding distance of a TM in the raw connectivity
    #   (e) Current bo: N_alpha-N_beta=1 (single) and N_beta-N_gamma=2 (double) — wrong order
    #
    # Fix: set bo(N_alpha, N_beta) = 2 (double), bo(N_beta, N_gamma) = 1 (single).
    # Then re-run upgrade_bonds to let the greedy add the SECOND double bond
    # (N_alpha=N_beta=N_gamma, or upgrade N_alpha-N_beta to triple for the other resonance).
    #
    # For CBC: the correct end-on azide N_alpha should appear as L (lone-pair donor).
    # With N_alpha=N_beta (bo=2): n_sub_organic=2, deficit=3-2=1 -> X (still wrong).
    # With N_alpha≡N_beta (bo=3): n_sub_organic=3, deficit=0, lp>0 -> L (CORRECT).
    # The second upgrade pass can achieve the triple if both N_alpha and N_beta are deficient.
    #
    _raw_adj_azide = defaultdict(set)
    for _a, _b in connectivity(atoms, coords):
        _raw_adj_azide[_a].add(_b); _raw_adj_azide[_b].add(_a)

    _azide_fixed = True
    while _azide_fixed:
        _azide_fixed = False
        for _nb in range(len(atoms)):          # N_beta candidate
            if atoms[_nb] != 'N':
                continue
            _nb_nbrs = list(adj[_nb])
            if len(_nb_nbrs) != 2:
                continue
            _na, _ng = _nb_nbrs[0], _nb_nbrs[1]
            # Ensure both ends are N
            if atoms[_na] != 'N' or atoms[_ng] != 'N':
                continue
            # Ensure both ends are terminal (exactly 1 Lewis bond each)
            if len(adj[_na]) != 1 or len(adj[_ng]) != 1:
                continue
            _key_ab = (min(_na, _nb), max(_na, _nb))
            _key_bg = (min(_nb, _ng), max(_nb, _ng))
            _bo_ab = bo.get(_key_ab, 1)
            _bo_bg = bo.get(_key_bg, 1)
            # Check if the wrong end has the double bond:
            # N_alpha (TM-adjacent, coordinating end) should have the HIGHER bond order.
            # Identify which terminal N is TM-adjacent.
            _na_has_tm = any(is_TM(atoms[_k]) for _k in _raw_adj_azide[_na])
            _ng_has_tm = any(is_TM(atoms[_k]) for _k in _raw_adj_azide[_ng])
            if _na_has_tm == _ng_has_tm:
                continue  # ambiguous or bridging azide — skip
            # _alpha is the TM-adjacent end, _gamma is the free terminal end
            _alpha, _gamma = (_na, _ng) if _na_has_tm else (_ng, _na)
            _key_alpha_beta = (min(_alpha, _nb), max(_alpha, _nb))
            _key_beta_gamma = (min(_nb, _gamma), max(_nb, _gamma))
            _bo_alpha = bo.get(_key_alpha_beta, 1)
            _bo_gamma = bo.get(_key_beta_gamma, 1)
            # Correct pattern: bo(alpha-beta) >= bo(beta-gamma)
            if _bo_alpha >= _bo_gamma:
                continue  # already correct or equal — nothing to swap
            # Wrong pattern: double bond is on the gamma side. Fix by swapping.
            bo[_key_alpha_beta] = _bo_gamma
            bo[_key_beta_gamma] = _bo_alpha
            _azide_fixed = True
            break  # restart scan

    # After swapping, re-run upgrade to let greedy add the second double / triple bond
    # on the alpha-beta pair (now with the single bond correctly on beta-gamma side).
    upgrade_bonds(atoms, coords, bo, adj)

    # Step 2: post-upgrade formal charge redistribution
    fc = redistribute_formal_charges(atoms, bo, adj, charge)

    # Step 3: lone pairs
    lp = compute_lone_pairs(atoms, bo, adj, fc)

    # Step 3a: Free-fragment lone-pair correction.
    #
    # $CHOOSE LONE must show the correct lone-pair count for each atom.
    # The key rule: for any atom, bond_e + 2*LP must equal VALENCE - fc.
    #
    # Two cases need correction:
    # A. Atom NOT sigma-bonded to TM in Lewis (dative-only or no contact):
    #    If total_e = bond_e + 2*lp < VALENCE[sym] AND fc=0 → assign fc=-1, recompute lp.
    #    Examples: N with 2 bonds (H+C), lp=1 computed → total=4 < 5 → fc=-1, lp=2 ✓
    #              Anionic O from aromatic ring enforcement with 1 bond, lp=2, total=5 < 6 ✓
    #
    # B. Atom sigma-bonded to TM in Lewis (bond_e_full includes TM sigma bond):
    #    If bond_e_full == VALENCE → atom is satisfied; keep fc=0, lp=0 as-is (correct).
    #    If bond_e_full < VALENCE → atom is STILL deficient even with TM bond → fc=-1.
    #    Note: for most X-type atoms (acyl C, alkynyl C, alkyl C bonded to TM),
    #    bond_e_full = VALENCE so they correctly have fc=0, lp=0.
    #    The rare case where bond_e_full < VALENCE (e.g. halide I with no Lewis bonds)
    #    is handled by fixing the ionic-metals path in redistribute_formal_charges.

    tm_set_3a = {i for i, sym in enumerate(atoms) if is_TM(sym)}

    for i, sym in enumerate(atoms):
        if is_TM(sym) or is_H(sym):
            continue
        if fc[i] != 0:
            continue   # fc already assigned — trust the redistribution result

        bond_e_full = sum(bo.get((min(i, k), max(i, k)), 0) for k in adj[i])
        bond_e_to_tm = sum(bo.get((min(i, k), max(i, k)), 0) for k in adj[i]
                           if k in tm_set_3a)

        exp_neutral = VALENCE.get(sym, 4)

        if bond_e_to_tm > 0:
            # Atom IS sigma-bonded to TM in Lewis.
            # Only correct if still deficient in the FULL Lewis structure.
            if bond_e_full < exp_neutral:
                fc[i] = -1
                eff_v = exp_neutral + 1
                lp[i] = max(0, eff_v - bond_e_full) // 2
            # If bond_e_full == exp_neutral: satisfied → keep fc=0, lp as-is.
        else:
            # Atom NOT sigma-bonded to TM (dative-only or isolated).
            total_e_now = bond_e_full + 2 * lp.get(i, 0)
            if total_e_now < exp_neutral:
                fc[i] = -1
                eff_v = exp_neutral + 1
                lp[i] = max(0, eff_v - bond_e_full) // 2


    # Step 3b: acyl / beta-diketonate LP shift
    # The greedy bond-upgrade sometimes places two C=O double bonds in a symmetric
    # acyl system (e.g. malonyl, beta-diketonate) leaving a lone pair on the central
    # carbon (C between two carbonyls) instead of on oxygen.
    # The more electronegative oxygen should carry the lone pair.
    # Fix: when a carbon has LP > 0, is NOT adjacent to any double bond of its own,
    # but IS adjacent to a carbonyl carbon (neighbour C has a C=O double bond),
    # shift one C=O → C-O single bond + add C=C double bond, transfer LP to O.
    _shift_acyl_lp(atoms, bo, adj, lp, fc)

    # ── Electron counts ───────────────────────────────────────────────────────
    core_e  = sum(CORE_E.get(a, 0) for a in atoms)
    val_e   = sum(VALENCE.get(a, 0) for a in atoms) - charge
    total_e = core_e + val_e

    n_sigma = len(bo)
    n_pi    = sum(o-1 for o in bo.values() if o >= 2)
    n_lp    = sum(lp.values())
    n_heavy = sum(1 for a in atoms if a != 'H')
    n_H     = atoms.count('H')

    val_nl  = NL_SIGMA*n_sigma + NL_PI*n_pi + NL_LP*n_lp
    ryd_nl  = RYD_HEAVY*n_heavy + RYD_H*n_H
    lewis_e = total_e - val_nl - ryd_nl

    return bo, lp, fc, dict(
        total_e=float(total_e), lewis_e=lewis_e,
        val_nl=val_nl, ryd_nl=ryd_nl, charge=float(charge)
    )


# ─── Output ───────────────────────────────────────────────────────────────────

def print_summary_and_choose(atoms, bo, lp, stats):
    te, le = stats['total_e'], stats['lewis_e']
    vnl, ryl = stats['val_nl'], stats['ryd_nl']
    print("          -------------------------------")
    print(f"                 Total Lewis{le:10.5f}  ({100*le/te:8.4f}%)")
    print(f"           Valence non-Lewis{vnl:10.5f}  ({100*vnl/te:8.4f}%)")
    print(f"           Rydberg non-Lewis{ryl:10.5f}  ({100*ryl/te:8.4f}%)")
    print("          -------------------------------")
    print(f"               Total unit  1{te:10.5f}  ({100.0:8.4f}%)")
    print(f"              Charge unit  1{stats['charge']:10.5f}")
    print()
    print(" $CHOOSE")
    lone_items = [(i+1, v) for i, v in sorted(lp.items()) if v > 0]
    if lone_items:
        print("   LONE " + " ".join(f"{idx} {v}" for idx,v in lone_items) + " END")
    sym = {1:'S', 2:'D', 3:'T'}
    tokens = [f"{sym[o]} {i+1} {j+1}" for (i,j), o in sorted(bo.items())]
    line = "   BOND"
    for tok in tokens:
        if len(line)+1+len(tok) > 72:
            print(line); line = "       "+tok
        else:
            line += " "+tok
    print(line + " END")
    print(" $END")


# ─── File readers ─────────────────────────────────────────────────────────────

def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    atoms, coords = [], []
    for ln in lines[2:2+n]:
        p = ln.split()
        raw = p[0]
        sym = raw[0].upper() + raw[1:].lower() if len(raw) > 1 else raw.upper()
        atoms.append(sym)
        coords.append([float(x) for x in p[1:4]])
    return atoms, coords, 0

_ANUM = {
     1:'H',  2:'He', 3:'Li', 4:'Be', 5:'B',  6:'C',  7:'N',  8:'O',
     9:'F', 10:'Ne',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P', 16:'S',
    17:'Cl',18:'Ar',19:'K', 20:'Ca',21:'Sc',22:'Ti',23:'V', 24:'Cr',
    25:'Mn',26:'Fe',27:'Co',28:'Ni',29:'Cu',30:'Zn',31:'Ga',32:'Ge',
    33:'As',34:'Se',35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y', 40:'Zr',
    41:'Nb',42:'Mo',43:'Tc',44:'Ru',45:'Rh',46:'Pd',47:'Ag',48:'Cd',
    49:'In',50:'Sn',51:'Sb',52:'Te',53:'I', 54:'Xe',55:'Cs',56:'Ba',
    57:'La',72:'Hf',73:'Ta',74:'W', 75:'Re',76:'Os',77:'Ir',78:'Pt',
    79:'Au',80:'Hg',81:'Tl',82:'Pb',83:'Bi',
}

def read_nbo47(path):
    with open(path) as f: text = f.read()
    atoms, coords, in_coord = [], [], False
    for line in text.splitlines():
        s = line.strip()
        if '$COORD' in s.upper(): in_coord = True; continue
        if in_coord and '$END' in s.upper(): break
        if in_coord and s and not s.startswith('$'):
            p = s.split()
            if len(p) >= 5:
                try:
                    sym = _ANUM.get(int(p[1]), f'X{p[1]}')
                    atoms.append(sym); coords.append([float(p[2]),float(p[3]),float(p[4])])
                except ValueError: pass
    return atoms, coords, 0


# ─── Demo: methylamine CH3NH2 (coords from nbo.log) ──────────────────────────
DEMO_ATOMS  = ['C','N','H','H','H','H','H']
DEMO_COORDS = [
    [ 0.745914,  0.011106,  0.000000],
    [-0.721743, -0.071848,  0.000000],
    [ 1.042059,  1.060105,  0.000000],
    [ 1.129298, -0.483355,  0.892539],
    [ 1.129298, -0.483355, -0.892539],
    [-1.076988,  0.386322, -0.827032],
    [-1.076988,  0.386322,  0.827032],
]

# ─── CBC Ligand Classification ────────────────────────────────────────────────
#
# The TRUE CBC rule (Green 1995, Green & Parkin 2014):
#
#   X-type : the ligand atom needs one more bond to complete its own normal
#             valence. It contributes ONE electron to the M-L bond (covalent).
#             Derived from ANIONIC precursors: H⁻, CH₃⁻, Cl⁻, OH⁻, NR₂⁻, OR⁻.
#             Test: bonds_to_substituents + 1(bond to M) == normal_valence → X.
#
#   L-type : the ligand atom already has its normal bond count satisfied by
#             its substituents ALONE. It contributes TWO electrons as a lone
#             pair (dative bond). Derived from NEUTRAL precursors: NH₃, NR₃,
#             OH₂, OR₂, CO, PR₃, alkenes, etc.
#             Test: bonds_to_substituents == normal_valence  AND  lone pair available.
#
#   Z-type : the ligand atom accepts two electrons FROM the metal (Lewis acid).
#             Rare. Examples: BR₃ with empty orbital, hypervalent I(III).
#
# For ambident ligands (NR₂ with lone pair AND needing one bond = amide):
#   → NR₂ with pyramidal N: X  (needs 1 bond, like all amides)
#   → NR₂ with planar N (π-donor): XL  (one bond needed + lone pair donation)
#   For simplicity we classify by bond-count alone (pyramidal approximation).
#
# Haptic π-systems:
#   alkene η²    → 1L
#   allyl⁻ η³   → LX
#   Cp⁻ η⁵      → L₂X
#   benzene η⁶  → L₃
#
# Atom categories (per user spec):
#   Organic-type : C, H, N, O, F, Cl  (appear as ligands to inorganic centres)
#   Inorganic-type: everything else — TMs, B, Si, P, S, Br, I, heavier halogens

_ORGANIC = {'C', 'H', 'N', 'O', 'F', 'Cl'}

def is_inorganic(sym):
    return sym not in _ORGANIC

# Normal (neutral) valence used for bond-count test
# For CBC we use the standard covalent bond number, not electron count
_NORMAL_BONDS = {
    'H':1,
    'B':3, 'C':4, 'N':3, 'O':2, 'F':1,
    'Si':4, 'P':3, 'S':2, 'Cl':1,
    'Ge':4, 'As':3, 'Se':2, 'Br':1,
    'Sn':4, 'Sb':3, 'Te':2, 'I':1,
    'Pb':4, 'Bi':3,
    'Al':3, 'Ga':3, 'In':3, 'Tl':3,
}

def normal_bonds(sym):
    """Standard number of covalent bonds in a neutral fragment."""
    return _NORMAL_BONDS.get(sym, 4)

def _subscript(n):
    subs = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
    return str(n).translate(subs) if n > 1 else ''


def classify_cbc_ligands(atoms, coords, bo, lp, fc, charge=0):
    """
    For each inorganic atom M, classify every bonded neighbour as L, X, or Z.

    Core logic for atom A bonded to metal M:
      n_sub = number of bonds A makes to atoms OTHER THAN M (from Lewis structure)
      n_val = normal_bonds(A)   ← bonds needed for neutral A to be complete

      if n_sub == n_val:          # already satisfied → lone-pair donor
          classify L  (dative)
      elif n_sub == n_val - 1:    # needs 1 more bond → covalent bond to M
          classify X
      elif n_sub <  n_val - 1:    # needs >1 bond (rare, e.g. bare atom) → X per bond needed
          classify X * (n_val - n_sub)   (each unsatisfied valence is one X)
      elif n_sub >  n_val:        # hypervalent substituent count
          → check lone pairs; if has LP → L; else Z

    Special cases handled separately:
      - π-system / haptic clusters (alkene, allyl, Cp, benzene)
      - CO / isocyanide (C is satisfied by triple bond to O but donates via LP → L)
      - Heavy halides Br, I: treated like Cl (X when simple halide; L when organic R-X)
      - TM neighbours: use electronegativity to decide L vs Z
      - Z-type: BR₃ (empty orbital), hypervalent I(III) with no LP residual
    """
    n = len(atoms)

    # ── Full connectivity for dative bond detection ──────────────────────────
    # For CBC classification we need to detect ALL TM–ligand contacts including:
    #   - Long dative TM–N/O bonds (up to 1.6× cov_r sum)
    #   - Haptic TM–C contacts in rings (up to 1.3× cov_r sum for ring members)
    #   - Standard sigma TM–C bonds (1.15×)
    # Strategy: use 1.15 for TM-C by default; then post-expand to include ring-C
    # atoms that are adjacent to already-bonded ring atoms (catch "missing" haptic C).
    def _full_connectivity_dative(atoms, coords, bo_lewis, lp_lewis):
        """
        Build the full connectivity list used for CBC dative-bond detection.

        Rules (applied in priority order):
        1. Standard 1.3× covalent-radii threshold for all bonds.
        2. For TM–heteroatom bonds (N, O, S, P, halogens) in the 1.3–1.6× zone:
           INCLUDE only if:
             a. The ligand atom has lp > 0 in the Lewis structure ($CHOOSE), AND
             b. The ligand atom is NOT bonded (in Lewis) to any atom that is
                already sigma-bonded to the TM (backbone exclusion: NHC N atoms,
                backbone amines, etc. that are 2 bonds away from TM via a sigma link
                are excluded — they are ligand framework, not additional donors), AND
             c. The ligand atom is NOT saturated with no lone pairs
                (bond_e >= std_cap AND lp == 0 → cannot donate).
        3. For TM–C bonds: 1.15× base; haptic ring C members expanded via BFS
           from dative C (non-σ) atoms.
        4. Terminal CO/CN oxygen exclusion (far end of already-coordinated ligand).
        5. Saturated atoms (bond_e >= std_cap AND lp == 0) are NEVER donors, even
           at standard 1.3× distance.
        """
        lewis_adj = defaultdict(set)
        for (a, b) in bo_lewis:
            lewis_adj[a].add(b); lewis_adj[b].add(a)

        # For each TM, collect atoms directly sigma-bonded to it in Lewis
        tm_sigma_partners = defaultdict(set)  # tm_idx → set of Lewis-bonded ligand atoms
        for a, b in bo_lewis:
            sa, sb = atoms[a], atoms[b]
            if is_TM(sa) and not is_TM(sb):
                tm_sigma_partners[a].add(b)
            elif is_TM(sb) and not is_TM(sa):
                tm_sigma_partners[b].add(a)

        def is_terminal_co_o(lig_idx, tm_idx):
            """Terminal O/N at the far end of a CO/CN ligand bonded via C to TM."""
            sym = atoms[lig_idx]
            if sym not in ('O', 'N'):
                return False
            nbrs = list(lewis_adj[lig_idx])
            if len(nbrs) != 1:
                return False
            c_nbr = nbrs[0]
            if atoms[c_nbr] != 'C':
                return False
            return tm_idx in lewis_adj[c_nbr]

        def is_backbone_atom(lig_idx, tm_idx):
            """
            Return True if lig_idx is a backbone heteroatom that should NOT be
            counted as an additional donor.
            A backbone atom is bonded (in Lewis) to at least one atom that is
            ITSELF directly sigma-bonded to the TM.
            Examples: NHC ring N atoms bonded to the carbene C which bonds to TM;
                      backbone amines two bonds away from coordinated donor.
            Only applies in the EXTENDED zone (1.3-1.6×). At ≤1.3× we trust geometry.
            """
            for nbr in lewis_adj[lig_idx]:
                if nbr in tm_sigma_partners[tm_idx]:
                    return True
            return False

        def is_agostic_h(lig_idx, tm_idx):
            """
            Return True if this H is part of an agostic or σ-complex interaction:
            H has exactly ONE Lewis bond to a C or B atom (not to the TM), and the
            H-to-TM distance is short enough to be a real agostic contact.
            These H atoms donate the B-H or C-H σ-bond to the metal → L-type.
            Also covers dihydrogen (H-H bond where both H are near the metal).
            """
            sym = atoms[lig_idx]
            if sym != 'H':
                return False
            nbrs = list(lewis_adj[lig_idx])
            if len(nbrs) == 0:
                # H with NO Lewis bond: could be dihydrogen partner — check if
                # there is another H very close (H-H dist < ~1.1 Å) and both near TM.
                ri_h = COV_R.get('H', 0.31)
                for other_h in range(len(atoms)):
                    if other_h == lig_idx or atoms[other_h] != 'H':
                        continue
                    d_hh = dist(coords[lig_idx], coords[other_h])
                    if d_hh < 1.15:   # H-H bond distance threshold
                        other_nbrs = list(lewis_adj[other_h])
                        if len(other_nbrs) == 0:   # other H also unbound in Lewis
                            return True   # dihydrogen η²-H₂
                return False
            if len(nbrs) == 1:
                nbr = nbrs[0]
                if is_TM(atoms[nbr]):
                    return False   # H bonded to TM in Lewis = classical hydride
                if atoms[nbr] in ('B', 'C', 'Si', 'Al', 'Ga'):
                    return True   # agostic B-H, C-H, Si-H → σ-complex
            return False

        def is_saturated_no_lp(lig_idx):
            """
            Return True if atom cannot coordinate: no lone pairs AND bond count
            at or beyond standard capacity.
            Carbon is NEVER excluded this way — it can be a haptic π-donor
            even with bond_e == 4 (the π-system donates, not a lone pair).
            H atoms that are part of agostic/σ-complex interactions are also
            included (handled separately by is_agostic_h).
            Only applies to heteroatoms (N, O, S, P, halogens).
            """
            sym = atoms[lig_idx]
            if sym == 'C':
                return False   # C can always be a haptic donor via π-system
            if is_TM(sym):
                return False
            if sym == 'H':
                # Agostic/dihydrogen H are never "saturated" for our purposes
                # — they are included via the agostic_h path.
                # Regular H (no Lewis bond or Lewis bond to TM) → exclude as always.
                return False   # let the TM-H distance check handle inclusion
            lp_count = lp_lewis.get(lig_idx, 0) if lp_lewis else 0
            if lp_count > 0:
                return False   # has lone pairs → can donate
            cap = STD_CAP.get(sym, 4)
            bond_e = sum(bo_lewis.get((min(lig_idx, k), max(lig_idx, k)), 0)
                         for k in lewis_adj[lig_idx])
            return bond_e >= cap   # fully saturated with no LP → cannot donate

        bonds = []
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                ri = COV_R.get(atoms[i], 0.77)
                rj = COV_R.get(atoms[j], 0.77)
                si, sj = atoms[i], atoms[j]
                tm_i, tm_j = is_TM(si), is_TM(sj)
                d_ij = dist(coords[i], coords[j])

                if (tm_i and not tm_j) or (tm_j and not tm_i):
                    tm = i if tm_i else j
                    lig = j if tm_i else i
                    sym_lig = atoms[lig]
                    lp_count = lp_lewis.get(lig, 0) if lp_lewis else 0

                    # Exclude terminal CO/CN oxygens
                    if is_terminal_co_o(lig, tm):
                        continue

                    # Exclude fully saturated atoms with no LP (cannot donate)
                    if is_saturated_no_lp(lig):
                        continue

                    threshold_std = 1.3 * (ri + rj)

                    # ── H atoms: agostic/σ-complex and dihydrogen special cases ──
                    if sym_lig == 'H':
                        # Only include H atoms that are agostic (bonded to B/C/Si/Al),
                        # dihydrogen (no Lewis bonds, or Lewis bond only to another H),
                        # or bare. Exclude protic H (bonded to O, N, S, F, etc.).
                        h_nbrs = list(lewis_adj[lig])
                        if len(h_nbrs) == 0:
                            # Bare H (no Lewis bond) — might be dihydrogen; include
                            if d_ij < threshold_std:
                                bonds.append((i, j))
                        elif len(h_nbrs) == 1:
                            parent_sym = atoms[h_nbrs[0]]
                            if parent_sym in ('B', 'C', 'Si', 'Al', 'Ga'):
                                # Agostic / σ-complex H — include
                                if d_ij < threshold_std:
                                    bonds.append((i, j))
                            elif parent_sym == 'H':
                                # Dihydrogen: H-H Lewis bond added by Step 1f.
                                # Include this H so the dihydrogen pair appears in adj_full.
                                if d_ij < threshold_std:
                                    bonds.append((i, j))
                            # else: protic H (bonded to O/N/S/F/etc.) → exclude
                        continue  # skip the rest of the heteroatom/C logic

                    if sym_lig in ('N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se', 'Te'):
                        if d_ij < threshold_std:
                            bonds.append((i, j))  # always include within standard
                        elif lp_count >= 1 and d_ij < 1.6 * (ri + rj):
                            # Extended zone checks:
                            # 1. Not a backbone atom (bonded to TM's sigma partner)
                            # 2. Absolute distance cap: real TM-L dative bonds rarely
                            #    exceed ~3.0A (even for heavy atoms).
                            #    Use 1.45× as the effective cap for the extended zone
                            #    — this handles Ni-O(carbonyl) at 2.66A (1.51×Ni+O)
                            #    while rejecting Pd-N at 2.95A when backbone check fails.
                            # 3. For Br/I with large cov_r, cap more tightly at 1.35×
                            if is_backbone_atom(lig, tm):
                                pass  # exclude
                            else:
                                sym_lig_upper = sym_lig
                                if sym_lig in ('Br', 'I'):
                                    # Heavy halides: max 1.35× (aryl halides far away)
                                    if d_ij < 1.35 * (ri + rj):
                                        bonds.append((i, j))
                                else:
                                    # N, O, S, P: max 1.55×
                                    if d_ij < 1.55 * (ri + rj):
                                        bonds.append((i, j))
                    elif sym_lig == 'C':
                        if d_ij < 1.15 * (ri + rj):
                            bonds.append((i, j))  # sigma TM-C
                        # haptic ring C added in post-expand below
                    else:
                        if d_ij < threshold_std:
                            bonds.append((i, j))
                else:
                    if d_ij < 1.3 * (ri + rj):
                        bonds.append((i, j))

        # Post-expand: for each TM, find its DATIVE C contacts (not σ-bonded in Lewis),
        # then BFS through C-C bonds to catch ring members at slightly longer distances.
        # Key: only expand from C atoms that are dative (no Lewis σ-bond to TM),
        # because σ-bonded C atoms (e.g. NHC carbenes, alkyls) are not haptic ring
        # members and their substituents must not be pulled into the TM neighbourhood.
        tm_indices = [i for i in range(len(atoms)) if is_TM(atoms[i])]
        bond_set = set(map(tuple, [sorted(b) for b in bonds]))
        adj_tmp = defaultdict(set)
        for i, j in bonds:
            adj_tmp[i].add(j); adj_tmp[j].add(i)

        # Build C-C adjacency for BFS
        all_C_adj = defaultdict(set)
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                if atoms[i] == 'C' and atoms[j] == 'C':
                    ri2 = COV_R.get(atoms[i], 0.77)
                    rj2 = COV_R.get(atoms[j], 0.77)
                    if dist(coords[i], coords[j]) < 1.35 * (ri2 + rj2):
                        all_C_adj[i].add(j); all_C_adj[j].add(i)

        for tm in tm_indices:
            # Only dative TM-C contacts are haptic seeds (no Lewis σ-bond)
            tm_C_all = {j for j in adj_tmp[tm] if atoms[j] == 'C'}
            tm_C_dative = {c for c in tm_C_all
                           if bo_lewis.get((min(tm, c), max(tm, c)), 0) == 0}
            if not tm_C_dative:
                continue

            # BFS from dative C atoms through C-C ring bonds
            visited_c = set(tm_C_all)   # avoid re-adding already known contacts
            frontier = set(tm_C_dative)
            for _ in range(3):  # up to 3 hops to cover 6-membered rings
                new_frontier = set()
                for c in frontier:
                    for nb in all_C_adj[c]:
                        if nb not in visited_c:
                            ri2 = COV_R.get(atoms[tm], 0.77)
                            rj2 = COV_R.get(atoms[nb], 0.77)
                            d_nb = dist(coords[tm], coords[nb])
                            if d_nb < 1.3 * (ri2 + rj2):
                                key = (min(tm, nb), max(tm, nb))
                                if key not in bond_set:
                                    bonds.append((tm, nb))
                                    bond_set.add(key)
                                    adj_tmp[tm].add(nb)
                                    adj_tmp[nb].add(tm)
                                visited_c.add(nb)
                                new_frontier.add(nb)
                frontier = new_frontier

        return bonds

    all_bonds_list = _full_connectivity_dative(atoms, coords, bo, lp)
    adj_full = defaultdict(set)
    for i, j in all_bonds_list:
        adj_full[i].add(j); adj_full[j].add(i)

    # ── Lewis connectivity (covalent only, dative bonds removed) ─────────────
    adj_lewis = defaultdict(set)
    for (i, j) in bo:
        adj_lewis[i].add(j); adj_lewis[j].add(i)

    def bo_ij(i, j):
        return bo.get((min(i, j), max(i, j)), 0)

    def n_sub_bonds(nbr_idx, metal_idx):
        """Total bond-order of nbr to all atoms EXCEPT the query metal."""
        return sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx] if k != metal_idx)

    # ── π-system detection ───────────────────────────────────────────────────

    def _pi_cluster(metal_idx, seed_idx):
        """
        Return the set of C or N atoms mutually connected via C–C or N–N or C–N bonds
        that are all bonded to metal_idx (haptic π-systems).

        Covers: alkene (C=C), alkyne, allyl, Cp, benzene (C-only clusters),
        AND: diazenido/azomethine/N=N ligands (N–N clusters),
             imine/pyridyl π-systems coordinating η² (C–N clusters).
        """
        if atoms[seed_idx] not in ('C', 'N'):
            return None
        # Metal's C/N neighbours
        metal_CN = {j for j in adj_full[metal_idx] if atoms[j] in ('C', 'N')}
        if seed_idx not in metal_CN:
            return None
        # BFS through C/N connections within the metal-bonded set
        visited = set()
        queue = [seed_idx]
        while queue:
            cur = queue.pop()
            if cur in visited: continue
            visited.add(cur)
            for nb in adj_full[cur]:
                if nb in metal_CN and atoms[nb] in ('C', 'N') and nb not in visited:
                    queue.append(nb)
        return visited

    def _classify_pi(metal_idx, nbr_idx):
        """
        If nbr_idx is part of a haptic π-system bound to metal_idx, return the
        CBC type list for the whole cluster (reported once for the lowest-index atom).
        Returns None if this is not a haptic case or not the representative atom.

        Single contacts (C or N, dative to M):
          - CO (C with C=O) → L
          - Alkene C (C=C, dative to M) → L
          - N with a π-bond (C=N, N=N, N=O) and dative to M → L (lone pair on N)
          - Everything else → None (handled by the bond-count rule)

        Multi-atom haptic clusters (representative = lowest index):
          C-only:  alkene η², allyl η³, Cp⁻ η⁵, benzene η⁶, etc.
          N-only or C/N mixed:  N=N η² (diazenido) → L, nitrile η², imine η², etc.
        """
        if atoms[nbr_idx] not in ('C', 'N'):
            return None

        cluster = _pi_cluster(metal_idx, nbr_idx)

        # ── Single contact ───────────────────────────────────────────────────
        if cluster is None or len(cluster) <= 1:
            # Only treat as π-donor when atom is NOT sigma-bonded to M in Lewis.
            if bo_ij(nbr_idx, metal_idx) > 0:
                return None  # covalent σ bond → handled by bond-count rule

            if atoms[nbr_idx] == 'C':
                # Guard: a C that already has its full complement of bonds to
                # non-metal substituents (bond_e >= normal_bonds) AND carries no
                # lone pairs cannot act as a lone-pair / π-donor.
                # This prevents fully-bonded ring carbons (e.g. indene/Cp junction
                # carbons with bond_e == 4) from being mis-classified as L donors
                # just because they happen to hold a C=C double bond.
                n_val_c = normal_bonds('C')   # 4
                bond_e_c = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx])
                lp_c = lp.get(nbr_idx, 0) if lp else 0
                if bond_e_c >= n_val_c and lp_c == 0:
                    return None  # fully saturated, no LP → cannot donate as solo L

                # CO/isocyanide: C with C=O → L
                if any(atoms[k] == 'O' and bo_ij(nbr_idx, k) >= 2
                       for k in adj_lewis[nbr_idx]):
                    return ['L']
                # Alkene C with π-bond to substituent → L
                if any(bo_ij(nbr_idx, k) > 1 for k in adj_lewis[nbr_idx]):
                    return ['L']
                return None  # ordinary C → bond-count rule

            if atoms[nbr_idx] == 'N':
                # N with a π-bond to substituent AND dative to M → L (lone-pair donor)
                # Covers: nitroso N (N=O), imine N (C=N), diazenido terminal N (N=N)
                #
                # EXCEPTION: anionic azide α-N (fc=-1, bonded only to another N via N=N).
                # The N=N double bond here is internal to the azide ligand, not donated
                # to the metal as a π interaction. The α-N coordinates as an X-type
                # anionic donor (like a halide), contributing one electron covalently.
                # Guard: if N has fc=-1 AND all its Lewis neighbours are N (no C, O, S),
                # it is a terminal azide/nitrene N → fall through to the bond-count rule.
                fc_n = fc[nbr_idx] if fc else 0
                if fc_n < 0:
                    all_nbrs_are_N = all(atoms[k] == 'N'
                                         for k in adj_lewis[nbr_idx])
                    if all_nbrs_are_N:
                        return None   # anionic azide α-N → bond-count rule → X
                if any(bo_ij(nbr_idx, k) > 1 for k in adj_lewis[nbr_idx]):
                    return ['L']
                return None  # ordinary N → bond-count rule

            return None

        # ── Multi-atom haptic cluster ────────────────────────────────────────
        # The representative is the lowest-index DATIVE atom in the cluster
        # (not sigma-bonded to the metal in the Lewis structure).
        # Sigma-bonded atoms (like a Cp carbon that also has a Lewis σ-bond)
        # are classified independently by the bond-count rule and must not
        # serve as cluster representatives — doing so causes duplicate rows.
        dative_members = {j for j in cluster if bo_ij(j, metal_idx) == 0}
        if not dative_members:
            return None  # all atoms sigma-bonded → no haptic cluster to report
        if nbr_idx != min(dative_members):
            return None   # not the representative dative atom → skip

        size = len(cluster)
        # Count atoms that need covalent bonding from M (X-type):
        # Either formally anionic (fc < 0) OR structurally deficient
        # (bond_e < expected valence, meaning Lewis couldn't satisfy them).
        # Both cases represent atoms that need one more electron → X type.
        def _is_x_type(j):
            """An atom is X-type if it is formally anionic, structurally deficient,
            OR has a covalent sigma bond to the metal (already counted in its bond_e
            but still contributes covalently to the M-L bond)."""
            fc_j = fc[j] if fc else 0
            if fc_j < 0:
                return True  # formally anionic
            if bo_ij(j, metal_idx) > 0:
                return True  # sigma-bonded to TM → X-type member of haptic cluster
            sym_j = atoms[j]
            exp_be = 4 if sym_j == 'C' else (3 if sym_j == 'N' else 2)
            bond_e_j = sum(bo_ij(j, k) for k in adj_lewis[j])
            return bond_e_j < exp_be  # structurally deficient

        n_X_atoms = sum(1 for j in cluster if _is_x_type(j))

        # Classify the π-system donation
        if size == 2:
            a2, b2 = sorted(cluster)
            pi_order = max(0, bo_ij(a2, b2) - 1)
            if pi_order == 0:
                # Single bond between the two atoms — NO π bond to donate as η².
                # Classify each atom independently (fall through, return None).
                return None
            elif n_X_atoms > 0:
                pi_types = ['L', 'X']     # anionic/deficient → η²-LX
            elif pi_order >= 2:
                pi_types = ['L', 'L']     # triple bond donating both π sets
            else:
                pi_types = ['L']          # double bond η² → 1 L

        elif size == 3:
            # Allyl η³: always anionic (allyl⁻) → LX, regardless of whether
            # the Lewis structure assigned fc=-1 to a ring atom.
            pi_types = ['L', 'X']

        elif size == 4:
            # Diene η⁴: neutral → L₂; LX+L when one deficient atom detected.
            pi_types = ['L', 'L', 'X'] if n_X_atoms > 0 else ['L', 'L']

        elif size == 5:
            # Cp η⁵ (cyclopentadienyl): always anionic (Cp⁻) → L₂X.
            # Even when the Lewis structure assigns fc=0 to all ring carbons
            # (because dative bonds were removed before formal-charge assignment),
            # the correct CBC result is L₂X — one carbon MUST be X-type.
            # Forcing one X here prevents a spurious solo-L from n_L_solo > 0
            # that would otherwise appear in _expand_cluster.
            pi_types = ['L', 'L', 'X']

        elif size == 6:
            pi_types = ['L', 'L', 'L']

        else:
            # General rule: odd-membered clusters are always anionic (one X);
            # even-membered clusters are neutral (all L pairs).
            # For odd sizes we guarantee at least one X to avoid spurious
            # solo-L entries caused by n_L_solo > 0.
            n_X_forced = size % 2          # 1 for odd, 0 for even
            n_X_actual = max(n_X_forced, n_X_atoms)
            n_L = (size - n_X_actual) // 2
            pi_types = ['L'] * n_L + ['X'] * n_X_actual

        # If the cluster representative has its own lone pair (e.g. N in N=N or N=C),
        # add an extra L for the lone-pair (σ) donation, separate from the π donation.
        # This covers: diazenido N59 with LP + N59=N60 π-bond both donating to metal.
        # Only add when: representative has lp > 0 AND bond to M is dative (not in Lewis).
        rep_lp = lp.get(nbr_idx, 0) if lp else 0
        rep_bond_to_M = bo_ij(nbr_idx, metal_idx)
        if rep_lp > 0 and rep_bond_to_M == 0 and 'X' not in pi_types:
            # Add one extra L for the lone-pair σ-donation (only when cluster is pure L)
            pi_types = ['L'] + pi_types

        return pi_types

    # ── Per-atom classifier ──────────────────────────────────────────────────

    def classify_one(metal_idx, nbr_idx):
        """Returns (types_list, cluster_set_or_None)."""
        sym_m = atoms[metal_idx]
        sym_n = atoms[nbr_idx]

        # ── H atom: hydride, agostic (B-H/C-H σ-complex), or dihydrogen ─────
        if sym_n == 'H':
            h_nbrs = list(adj_lewis[nbr_idx])   # Lewis bonds of this H
            # Case 1: H bonded to TM in Lewis → classical hydride → X
            if any(is_TM(atoms[k]) for k in h_nbrs):
                return ['X'], None

            # Case 2: H bonded to another H (dihydrogen η²-H₂ ligand) AND
            # both close to the metal → L  (Only report for the lower-index H.)
            # After Step 1f in find_best_lewis the H-H σ-bond IS in the Lewis
            # structure, so h_nbrs may be [partner_H] rather than empty.
            # We detect dihydrogen by: H has exactly one Lewis bond and that
            # bond is to another H (not TM, not heavy atom).
            is_dihydrogen_h = False
            dihydrogen_partner = None
            if len(h_nbrs) == 0:
                # No Lewis bonds yet — search for a partner H within 1.15 Å
                for other_h in sorted(adj_full[metal_idx]):
                    if atoms[other_h] != 'H' or other_h == nbr_idx:
                        continue
                    other_nbrs = list(adj_lewis[other_h])
                    if any(is_TM(atoms[k]) for k in other_nbrs):
                        continue   # other H is a classical hydride
                    # Allow other_h to have at most one Lewis bond (to nbr_idx itself)
                    non_tm_nbrs = [k for k in other_nbrs if not is_TM(atoms[k])]
                    if len(non_tm_nbrs) > 1:
                        continue   # bonded to heavy atom → not dihydrogen
                    if len(non_tm_nbrs) == 1 and non_tm_nbrs[0] != nbr_idx:
                        continue   # bonded to a different atom → not dihydrogen
                    d_hh = dist(coords[nbr_idx], coords[other_h])
                    if d_hh < 1.15:
                        is_dihydrogen_h = True
                        dihydrogen_partner = other_h
                        break
            elif (len(h_nbrs) == 1 and atoms[h_nbrs[0]] == 'H'
                  and not any(is_TM(atoms[k]) for k in h_nbrs)):
                # H has a single Lewis bond to another H → dihydrogen (Step 1f added it)
                other_h = h_nbrs[0]
                if other_h in adj_full[metal_idx]:
                    is_dihydrogen_h = True
                    dihydrogen_partner = other_h

            if is_dihydrogen_h:
                if nbr_idx < dihydrogen_partner:
                    return ['L'], None   # lower-index H represents the pair
                else:
                    return [], None      # higher-index H: skip, already reported

            # Case 3: H bonded to B, C, Si, Al (agostic / σ-complex) → L
            # The B-H or C-H σ-bond donates to the metal as a dative interaction.
            sigma_donor_elements = ('B', 'C', 'Si', 'Al', 'Ga')
            if len(h_nbrs) == 1 and atoms[h_nbrs[0]] in sigma_donor_elements:
                # Confirm the bonded atom (B/C) is not itself sigma-bonded to this
                # metal in Lewis (that would make it a normal covalent ligand, not
                # agostic). If the parent is in adj_lewis[metal_idx] it's a standard
                # X-type and the H is just a substituent — exclude it.
                parent = h_nbrs[0]
                if bo_ij(parent, metal_idx) > 0:
                    # Parent is covalently bound to metal (in Lewis).
                    # If it's a C-H agostic (C sigma-bonded to M AND H near M),
                    # the interaction is agostic: the whole C-H σ-bond donates.
                    # In CBC this is classified as a separate L interaction.
                    # We include H as L only when H is actually within the standard
                    # dative distance of the metal (already confirmed by adj_full).
                    return ['L'], None
                else:
                    # Parent NOT in Lewis σ-bond to M: pure σ-complex (B-H, C-H)
                    return ['L'], None

            # Case 4: H bonded to N/O/S (protic H) near metal → X
            return ['X'], None

        # ── TM–TM bond ────────────────────────────────────────────────────────
        if is_TM(sym_n):
            en_m = _OXSTATE_ENEG.get(sym_m, 2.0)
            en_n = _OXSTATE_ENEG.get(sym_n, 2.0)
            return (['L'] if en_n < en_m else ['Z']), None

        # ── Inorganic main-group neighbour (B, Si, P, S, Br, I, etc.) ────────
        if is_inorganic(sym_n):
            n_val    = normal_bonds(sym_n)
            n_sub    = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx] if k != metal_idx)
            lp_lewis = lp.get(nbr_idx, 0)
            deficit  = n_val - n_sub

            if sym_n in ('Br', 'I') and n_sub == 0:
                return ['X'], None

            if sym_n == 'I' and n_sub >= 2 and is_TM(sym_m):
                fc_n = fc[nbr_idx] if fc else 0
                if fc_n > 0:
                    return ['Z'], None

            if sym_n in ('B', 'Al', 'Ga', 'In') and lp_lewis == 0:
                return ['Z'], None

            if deficit <= 0:
                return (['L'] if lp_lewis > 0 else ['Z']), None
            elif deficit == 1:
                return ['X'], None
            else:
                if lp_lewis > 0:
                    return (['L'] + ['X'] * (deficit - 1)), None
                else:
                    return (['X'] * deficit), None

        # ── Organic-type neighbour (C, H, N, O, F, Cl) ───────────────────────

        # First check for haptic π-system (only for C/N atoms)
        pi = _classify_pi(metal_idx, nbr_idx)
        if pi is not None:
            # Determine cluster for proper pair-row display.
            # ONLY attach cluster when _classify_pi confirmed haptic donation.
            cluster = _pi_cluster(metal_idx, nbr_idx)
            return pi, cluster

        # If this C or N is part of a multi-atom haptic cluster but is NOT the
        # lowest-index DATIVE representative, skip it entirely.
        # The representative is the lowest-index atom with no Lewis σ-bond to M.
        if atoms[nbr_idx] in ('C', 'N'):
            cluster_check = _pi_cluster(metal_idx, nbr_idx)
            if cluster_check and len(cluster_check) > 1:
                dative_check = {j for j in cluster_check if bo_ij(j, metal_idx) == 0}
                if dative_check and nbr_idx != min(dative_check):
                    return [], None  # non-representative dative haptic member — skip

        sym_n    = atoms[nbr_idx]
        n_val    = normal_bonds(sym_n)
        lp_lewis = lp.get(nbr_idx, 0)

        # n_sub_organic = bond-order sum to ORGANIC-TYPE neighbours only (excluding M
        # and inorganic heteroatoms like S, P, Si — dative bonds to other centres).
        n_sub_organic = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx]
                            if not is_inorganic(atoms[k]))

        # n_sub_all = bond-order sum to ALL non-metal neighbours (includes halides,
        # inorganic substituents that form real covalent bonds to C like C-I, C-Cl).
        # Used for carbene detection: a C bonded to I, H, H has 3 real bonds → NOT carbene.
        n_sub_all = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx]
                        if k != metal_idx and not is_TM(atoms[k]))

        deficit = n_val - n_sub_organic
        bond_in_lewis = bo_ij(nbr_idx, metal_idx) > 0

        # ── Special detection: CO and isocyanide ─────────────────────────────
        # CO is ONLY a lone-pair σ-donor via C when C has NO other organic bonds.
        # A C bonded to O (via C=O or C≡O) AND to other carbon-chain atoms is
        # NOT a terminal CO — it's a carbonyl in a larger molecule → use bond-count rule.
        if sym_n == 'C':
            has_CO = any(atoms[k] == 'O' and bo_ij(nbr_idx, k) >= 2
                         for k in adj_lewis[nbr_idx])
            if has_CO:
                # Only classify as CO if no other organic substituents
                other_organic = sum(1 for k in adj_lewis[nbr_idx]
                                    if k != metal_idx and not is_inorganic(atoms[k])
                                    and atoms[k] != 'O')
                if other_organic == 0:
                    return ['L'], None   # true terminal CO/isocyanide → L

        # ── Special detection: carbene ────────────────────────────────────────
        # Carbene: C with ≤ 2 bonds to ALL non-metal neighbours (not just organic).
        # A C-I or C-Cl bond is a real covalent bond that satisfies C's valence;
        # n_sub_all correctly counts it, preventing misclassification as carbene.
        if sym_n == 'C':
            if n_sub_all <= 2:
                return ['L'], None   # true carbene: ≤ 2 non-metal bonds → L

        # ── Main bond-count decision ──────────────────────────────────────────
        if deficit <= 0:
            if lp_lewis > 0:
                return ['L'], None
            return [], None          # saturated, no LP — cannot bind

        elif deficit == 1:
            return ['X'], None       # one organic bond short → covalent σ to M

        else:
            # deficit >= 2: atom has fewer organic bonds than normal valence.
            #
            # Sub-case (a): TERMINAL nitrido/oxo/imido atoms that are genuinely bare
            # (no bonds to any non-TM atom, or their only bond is to this metal and they
            # still need more bonds). The metal supplies all bonds covalently.
            # Each unsatisfied valence counts as one X-type interaction:
            #   nitrido N  (N≡M): n_sub_all_total=0, deficit_real=3 → 3X
            #   oxo     O  (O=M): n_sub_all_total=0, deficit_real=2 → 2X
            #   imido   NR with no further substituents: n_sub_all=0 except M bond
            #
            # IMPORTANT: we use n_sub_all_total (bonds to ALL non-TM atoms, INCLUDING
            # the Lewis bond to the metal center itself) to compute the real deficit.
            # This correctly distinguishes:
            #   - True nitrido/oxo: n_sub_all_total=0 (or just the M bond) → multi-X
            #   - Sulfonimide N (S-N-S): n_sub_all_total=2 → deficit_real=1 → 1X
            #   - Sulfonate O (S-O, only bond): n_sub_all_total=1 → deficit_real=1 → 1X
            # Using only n_sub_organic (which excludes inorganic S, P, Si) would
            # wrongly classify these as multi-X since S is inorganic.
            #
            # Sub-case (b): partially substituted atoms with real bonds to other
            # non-metals: use the real deficit for the multi-X test.
            n_sub_all_total = sum(bo_ij(nbr_idx, k) for k in adj_lewis[nbr_idx]
                                  if not is_TM(atoms[k]))
            deficit_real = n_val - n_sub_all_total
            if n_sub_all_total == 0 and sym_n in ('N', 'O', 'S', 'C'):
                # Truly terminal with no bonds to any non-TM atom at all —
                # all valence bonds go to the metal covalently
                return ['X'] * deficit, None
            elif deficit_real > 1 and sym_n in ('N', 'O', 'S', 'C'):
                # Still needs more than 1 bond even counting ALL its bonds
                # (e.g. imido NR with only 1 substituent, still 1 short for a true double):
                # use the real deficit, but cap at the organic deficit
                return ['X'] * min(deficit_real, deficit), None
            else:
                return ['X'], None


    def _fix_symmetric_chelate_O(metal_idx, nbr_cls, atoms, bo, adj_lewis, adj_full):
        """
        Post-processing correction for symmetric bidentate O chelates AND
        monodentate/bidentate oxyanion ligands (nitrate, carbonate, sulfate, etc.).

        Problem 1 — symmetric carboxylate/diketonate:
          Symmetric Lewis structure (both C=O double bonds) makes both O appear L.
          Correct CBC: one O is X (anionic, single bond to C), one O is L.

        Problem 2 — nitrate κ1 (monodentate) / κ2 (bidentate):
          One O bonds to M. That O has 1 bond to N → deficit=1 → X.  But N carries
          a double bond to another O in the ring or in the resonance structure.
          That second O (N=O) can donate a lone pair → L.
          A third O in κ1 is not coordinated at all.
          Correct CBC for κ1-NO3: 1X (coordinating O) + 1L (N=O oxygen lone pair
          donated through the N=O π system, or via N lone pair).

          More precisely: for oxyanion ligands (O atoms bonded to the same N or S or
          Cl central atom), when one O coordinates as X and another O of the same
          oxyanion is classified as L, the pair represents the correct κ2 coordination.
          For κ1, only one O of the pair actually coordinates; the L entry reflects
          the N=O resonance donation.

        Additionally fix symmetric carboxylate chelate pairs as before.
        """
        # ── 1. Carboxylate / diketonate symmetric fix (both O classified L via C=O) ──
        o_L = [idx for idx, types in nbr_cls.items()
               if atoms[idx] == 'O'
               and types == ['L']
               and bo.get((min(idx, metal_idx), max(idx, metal_idx)), 0) == 0]

        if len(o_L) >= 2:
            for i in range(len(o_L)):
                for j in range(i + 1, len(o_L)):
                    oi, oj = o_L[i], o_L[j]
                    # Both must have a C=O (bo≥2) to their carbon in Lewis
                    ci = next((k for k in adj_lewis[oi] if atoms[k] == 'C'
                               and bo.get((min(oi, k), max(oi, k)), 0) >= 2), None)
                    cj = next((k for k in adj_lewis[oj] if atoms[k] == 'C'
                               and bo.get((min(oj, k), max(oj, k)), 0) >= 2), None)
                    if ci is None or cj is None:
                        continue
                    connected = False
                    if ci == cj:
                        connected = True
                    elif cj in adj_lewis[ci] or ci in adj_lewis[cj]:
                        connected = True
                    else:
                        for mid in adj_lewis[ci]:
                            if cj in adj_lewis[mid]:
                                connected = True
                                break
                    if connected:
                        nbr_cls[max(oi, oj)] = ['X']
                        break

        # ── 3. Bidentate oxyanion fix (κ²-NO₃, κ²-CO₃, κ²-SO₄, etc.) ──────────────
        # When two or more O atoms of the same oxyanion (bonded to the same N/S/P/C
        # central atom) both coordinate to the metal and are both classified X, the
        # correct CBC treatment is: one O is X (anionic), one is L (dative lone-pair
        # donor via resonance). Reclassify the higher-index one as L.
        # Detection: two X-type O neighbours bonded to the same non-TM, non-metal
        # central atom (N, S, P, C) that itself is NOT directly bonded to the metal.
        x_o_nbrs = [idx for idx, types in nbr_cls.items()
                    if atoms[idx] == 'O' and types == ['X']
                    and bo.get((min(idx, metal_idx), max(idx, metal_idx)), 0) == 0]

        # Group these O atoms by their central heteroatom
        central_to_os = defaultdict(list)
        for o_idx in x_o_nbrs:
            for central_k in adj_lewis[o_idx]:
                if atoms[central_k] in ('N', 'S', 'P', 'C', 'Cl', 'Br') and \
                        not is_TM(atoms[central_k]) and central_k not in adj_full[metal_idx]:
                    central_to_os[central_k].append(o_idx)

        for central_k, o_list in central_to_os.items():
            if len(o_list) >= 2:
                # Reclassify the higher-index O(s) as L (keep the lowest as X)
                o_list_sorted = sorted(o_list)
                for o_to_flip in o_list_sorted[1:]:
                    # Only flip if this O has lone pairs (can donate)
                    if lp.get(o_to_flip, 0) > 0:
                        nbr_cls[o_to_flip] = ['L']
        # For each O neighbour classified as X (deficit=1, bonds to N/S/Cl/P heteroatom
        # via a single N-O bond), check whether another O of the SAME oxyanion group is
        # also in the metal neighbourhood and has a double bond to the central atom.
        # If so, that double-bonded O can donate a lone pair → classify it as L.
        # This correctly models κ1-NO3 as 1X + 1L.
        for x_idx, x_types in list(nbr_cls.items()):
            if atoms[x_idx] != 'O' or x_types != ['X']:
                continue
            # Find the central heteroatom this O is bonded to in Lewis (must be N, S, P, Cl)
            central = next((k for k in adj_lewis[x_idx]
                            if atoms[k] in ('N', 'S', 'P', 'Cl', 'Br', 'C')
                            and not is_TM(atoms[k])), None)
            if central is None:
                continue
            # The central atom must itself be bonded to ≥1 other O with a double bond
            for other_o in adj_lewis[central]:
                if other_o == x_idx or atoms[other_o] != 'O':
                    continue
                key_co = (min(central, other_o), max(central, other_o))
                if bo.get(key_co, 0) < 2:
                    continue
                # other_o has a double bond to the central atom → potential L donor
                # It must NOT already be classified (not in nbr_cls means it's not a
                # direct metal neighbour, which is the κ1 case; if it IS a neighbour
                # it would be in nbr_cls already and we just need to ensure it's L).
                if other_o in nbr_cls:
                    # Already classified: ensure it's L (double-bonded O is lone-pair donor)
                    if nbr_cls[other_o] == ['X'] and bo.get(key_co, 0) >= 2:
                        nbr_cls[other_o] = ['L']
                        break

    # ── Helper: expand haptic cluster types into interaction records ─────────
    def _expand_cluster(rep_idx, cluster_types, cluster_atoms_sorted, metal_idx):
        """
        Convert a flat list of CBC types (e.g. ['L','L','X']) for a haptic cluster
        into a list of interaction records using π-bond-guided pairing:
          - Pairs are formed from atoms sharing a double bond in the Lewis structure
            ($CHOOSE), ensuring chemically meaningful η² pairs.
          - Anionic atoms (fc < 0) are reserved for X rows.
          - A leading solo-atom L (lone-pair donation from rep atom) uses just rep_idx.
        """
        records = []
        types_work = list(cluster_types)
        n_atoms = len(cluster_atoms_sorted)
        n_L = types_work.count('L')
        n_L_cluster = n_atoms // 2
        n_L_solo = n_L - n_L_cluster

        # Emit solo LP records first (lone-pair σ-donation from rep atom)
        for _ in range(n_L_solo):
            records.append(((rep_idx,), 'L'))

        cluster_set = set(cluster_atoms_sorted)
        # Identify X-type atoms: formally anionic (fc<0) OR structurally deficient
        def _is_x_candidate(idx):
            fc_val = fc[idx] if fc else 0
            if fc_val < 0:
                return True
            sym = atoms[idx]
            exp_be = 4 if sym == 'C' else (3 if sym == 'N' else 2)
            bond_e = sum(bo_ij(idx, k) for k in adj_lewis[idx])
            return bond_e < exp_be

        anionic = {a for a in cluster_set if _is_x_candidate(a)}
        neutral = [a for a in cluster_atoms_sorted if a not in anionic]

        # Build π-bond pairs: find double bonds within cluster neutral atoms
        double_pairs = []
        used = set()
        # Collect all intra-cluster double bonds among neutral atoms
        bond_candidates = []
        for a in neutral:
            for b in adj_lewis[a]:
                if b in cluster_set and b not in anionic and b > a:
                    order = bo_ij(a, b)
                    bond_candidates.append((-order, a, b))   # higher order first
        bond_candidates.sort()
        for _, a, b in bond_candidates:
            if a not in used and b not in used and len(double_pairs) < n_L_cluster:
                double_pairs.append((a, b))
                used.add(a); used.add(b)

        # If not enough pairs from double bonds, use any adjacent neutral pair
        remaining = [a for a in neutral if a not in used]
        while len(double_pairs) < n_L_cluster and len(remaining) >= 2:
            a = remaining.pop(0)
            # Prefer an adjacent partner
            partner = next((b for b in remaining if b in adj_lewis[a]), None)
            if partner is None:
                partner = remaining[0]
            double_pairs.append((a, partner))
            remaining.remove(partner)
            used.add(a); used.add(partner)

        # Emit L pair rows
        for a, b in double_pairs:
            records.append(((a, b), 'L'))

        # Emit X rows: anionic atoms first, then any leftover neutrals
        x_atoms = sorted(anionic) + [a for a in cluster_atoms_sorted
                                      if a not in used and a not in anionic]
        for ax in x_atoms:
            records.append(((ax,), 'X'))

        return records


    # ── Main loop: report for each inorganic centre ──────────────────────────
    results = {}
    for metal_idx in range(n):
        sym_m = atoms[metal_idx]
        if not is_inorganic(sym_m):
            continue

        # Skip simple halides (Br, I with no own substituents) as stand-alone centres
        if sym_m in ('Br', 'I'):
            n_sub_m = sum(bo_ij(metal_idx, k) for k in adj_lewis[metal_idx])
            tm_nbrs = [j for j in adj_full[metal_idx] if is_TM(atoms[j])]
            if n_sub_m == 0 and tm_nbrs:
                continue

        neighbours = sorted(adj_full[metal_idx])
        if not neighbours:
            continue

        nbr_cls = {}       # nbr_idx -> simple ['L'],['X'],['Z'] for non-haptic
        nbr_cluster = {}   # nbr_idx -> cluster set (for haptic pi-donors)
        for nbr_idx in neighbours:
            cbc, cluster = classify_one(metal_idx, nbr_idx)
            if cbc:
                nbr_cls[nbr_idx] = cbc
                if cluster and len(cluster) > 1:
                    nbr_cluster[nbr_idx] = cluster

        _fix_symmetric_chelate_O(metal_idx, nbr_cls, atoms, bo, adj_lewis, adj_full)

        if not nbr_cls:
            continue

        # ── Expand haptic clusters into interaction records ───────────────────
        # interaction_records: list of (atom_tuple, cbc_char) ordered for display
        interaction_records = []
        reported_haptic_reps = set()

        for nbr_idx in sorted(nbr_cls.keys()):
            types = nbr_cls[nbr_idx]
            cluster = nbr_cluster.get(nbr_idx)

            # If this nbr has a multi-atom cluster (haptic pi-donor), expand to pair rows
            if cluster and len(cluster) > 1:
                cluster_sorted = sorted(cluster)
                records = _expand_cluster(nbr_idx, types, cluster_sorted, metal_idx)
                interaction_records.extend(records)
                reported_haptic_reps.add(nbr_idx)
                continue

            # ── Dihydrogen: show as a pair row H_i-H_j ────────────────────────
            if atoms[nbr_idx] == 'H' and types == ['L']:
                # Check if this H is the lower-index of a dihydrogen pair.
                # After Step 1f, dihydrogen H atoms have a mutual Lewis bond (H-H),
                # so h_nbrs_lewis may be [partner_H] instead of empty.
                h_nbrs_lewis = list(adj_lewis[nbr_idx])
                # Dihydrogen condition: H has no Lewis bonds OR its only Lewis bond
                # is to another H (not a heavy atom or TM).
                is_dh_candidate = (len(h_nbrs_lewis) == 0 or
                                   (len(h_nbrs_lewis) == 1
                                    and atoms[h_nbrs_lewis[0]] == 'H'
                                    and not any(is_TM(atoms[k]) for k in h_nbrs_lewis)))
                if is_dh_candidate:
                    # Find H-H partner in the metal neighbourhood
                    partner_h = None
                    for other_h in sorted(adj_full[metal_idx]):
                        if other_h <= nbr_idx or atoms[other_h] != 'H':
                            continue
                        other_nbrs_lewis = list(adj_lewis[other_h])
                        # Partner must also be dihydrogen-like (no Lewis bonds or
                        # only a Lewis bond to nbr_idx itself)
                        other_non_tm = [k for k in other_nbrs_lewis
                                        if not is_TM(atoms[k])]
                        if len(other_non_tm) > 1:
                            continue
                        if len(other_non_tm) == 1 and other_non_tm[0] != nbr_idx:
                            continue   # bonded to a different heavy atom
                        if any(is_TM(atoms[k]) for k in other_nbrs_lewis):
                            continue   # classical hydride
                        d_hh = dist(coords[nbr_idx], coords[other_h])
                        if d_hh < 1.15:
                            partner_h = other_h
                            break
                    if partner_h is not None:
                        # Emit as pair (H_i, H_j) L — dihydrogen
                        interaction_records.append(((nbr_idx, partner_h), 'L'))
                        continue
                # Agostic B-H / C-H / Si-H σ-complex: emit as (parent, H) pair
                # so the label reads "B2-H3" or "C6-H30" — the donation originates
                # from the σ-bond midpoint, not from the H atom alone.
                h_lewis_nbrs_agostic = list(adj_lewis[nbr_idx])
                if (len(h_lewis_nbrs_agostic) == 1
                        and atoms[h_lewis_nbrs_agostic[0]] in ('B', 'C', 'Si', 'Al', 'Ga')):
                    parent_idx = h_lewis_nbrs_agostic[0]
                    interaction_records.append(((parent_idx, nbr_idx), 'L'))
                else:
                    interaction_records.append(((nbr_idx,), types[0]))
                continue

            if len(types) == 1:
                # Simple single-type: one row
                interaction_records.append(((nbr_idx,), types[0]))
            elif len(set(types)) == 1 and len(types) > 1:
                # Multiple identical types (nitrido XXX, oxo XX): one row per bond
                for t in types:
                    interaction_records.append(((nbr_idx,), t))
            else:
                # Multi-type: expand into pair-rows (legacy path for non-cluster multi-type)
                cluster2 = _pi_cluster(metal_idx, nbr_idx)
                if cluster2 and len(cluster2) > 1 and nbr_idx == min(cluster2):
                    cluster_sorted = sorted(cluster2)
                    records = _expand_cluster(nbr_idx, types, cluster_sorted, metal_idx)
                    interaction_records.extend(records)
                    reported_haptic_reps.add(nbr_idx)
                else:
                    # Fallback for non-cluster multi-type (carbene LX, etc.)
                    for t in types:
                        interaction_records.append(((nbr_idx,), t))

        results[metal_idx] = interaction_records

    return results


# Electronegativity table for TM oxidation-state comparisons
_OXSTATE_ENEG = dict(ENEG)
_OXSTATE_ENEG.update({
    'Sc':1.36,'Ti':1.54,'V':1.63,'Cr':1.66,'Mn':1.55,'Fe':1.83,'Co':1.88,
    'Ni':1.91,'Cu':1.90,'Zn':1.65,'Y':1.22,'Zr':1.33,'Nb':1.60,'Mo':2.16,
    'Tc':1.90,'Ru':2.20,'Rh':2.28,'Pd':2.20,'Ag':1.93,'Cd':1.69,
    'Hf':1.30,'Ta':1.50,'W':2.36,'Re':1.90,'Os':2.20,'Ir':2.20,
    'Pt':2.28,'Au':2.54,'Hg':2.00,
})


def print_cbc_report(atoms, coords, bo, lp, fc, charge=0):
    """Print the CBC classification table for every inorganic centre.

    Results are shown as interaction records, one row per bond/donation:
      Single atom:   N1        2.05Å   L    N donates lone pair
      Pair (π bond): N59-N60   1.94Å   L    η²-N=N π donation
      Haptic ring:   C1-C2     2.5Å    L    one π bond of Cp ring
                     C3-C4     2.5Å    L
                     C5        2.5Å    X    anionic carbon
    """
    results = classify_cbc_ligands(atoms, coords, bo, lp, fc, charge)

    if not results:
        print("\n  (No inorganic atoms found — no CBC classification.)")
        return

    print()
    print("=" * 72)
    print("  CBC Ligand Classification  (L = dative, X = covalent, Z = Lewis-acid)")
    print("=" * 72)

    adj_lewis = defaultdict(set)
    for (i, j) in bo:
        adj_lewis[i].add(j); adj_lewis[j].add(i)

    def bo_ij(i, j):
        return bo.get((min(i, j), max(i, j)), 0)

    def avg_dist(metal_idx, atom_tuple):
        """Distance from metal to the interaction point.

        For single-atom rows: distance to that atom.
        For H-H dihydrogen pairs: average (midpoint of H-H to metal).
        For agostic X-H pairs (parent + H): distance to the H atom, which is
            the interaction point closest to the metal and the physically
            meaningful distance for the σ-complex contact.
        For π-bond pairs (C-C, N-N, etc.): average (midpoint of bond to metal).
        """
        if len(atom_tuple) == 1:
            a = atom_tuple[0]
            return math.sqrt(sum((coords[metal_idx][k]-coords[a][k])**2
                                 for k in range(3)))
        a, b = atom_tuple
        # Agostic X-H pair: one is H, one is B/C/Si/Al/Ga → use H distance
        if atoms[b] == 'H' and atoms[a] in ('B', 'C', 'Si', 'Al', 'Ga'):
            return math.sqrt(sum((coords[metal_idx][k]-coords[b][k])**2
                                 for k in range(3)))
        # All other pairs (H-H dihydrogen, C-C π-bond, etc.): midpoint
        dists = [math.sqrt(sum((coords[metal_idx][k]-coords[x][k])**2
                               for k in range(3))) for x in atom_tuple]
        return sum(dists) / len(dists)

    def row_label(atom_tuple):
        """Build the neighbour label: 'N59' or 'N59-N60' or 'C12-C14'."""
        parts = [f"{atoms[a]}{a+1}" for a in atom_tuple]
        return "-".join(parts)

    def row_note(atom_tuple, cbc_char, metal_idx):
        """Short explanatory note for one interaction row."""
        if cbc_char == 'Z':
            return "Lewis acid — metal donates electrons to ligand"
        rep = atom_tuple[0]
        sym = atoms[rep]
        n_val = _NORMAL_BONDS.get(sym, 4)
        n_sub = sum(bo_ij(rep, k) for k in adj_lewis[rep]
                    if not is_inorganic(atoms[k]))
        lp_cnt = lp.get(rep, 0) if lp else 0

        if len(atom_tuple) == 2:
            a, b = atom_tuple
            # Special case: H-H dihydrogen pair
            if atoms[a] == 'H' and atoms[b] == 'H':
                return "η²-H₂ dihydrogen — σ(H-H) dative donor"
            # Special case: agostic X-H σ-complex pair (parent atom + H)
            if atoms[b] == 'H' and atoms[a] in ('B', 'C', 'Si', 'Al', 'Ga'):
                return f"agostic {atoms[a]}-H σ-bond → σ-complex dative donation to M"
            # Multi-atom (π-bond pair): always a dative pi-bond donation
            bond_type = "C=C" if (atoms[a]=='C' and atoms[b]=='C') else (
                        "N=N" if (atoms[a]=='N' and atoms[b]=='N') else "π")
            return f"η² {bond_type} π-bond → dative π donation to M"

        # Single-atom interaction
        if sym == 'H':
            h_lewis_nbrs = list(adj_lewis[rep])
            if cbc_char == 'X':
                if h_lewis_nbrs:
                    parent = h_lewis_nbrs[0]
                    return f"hydride H — covalent M-H bond"
                return "hydride H — covalent M-H bond"
            if cbc_char == 'L':
                h_lewis_nbrs = list(adj_lewis[rep])
                if len(h_lewis_nbrs) == 1:
                    parent = h_lewis_nbrs[0]
                    parent_sym = atoms[parent]
                    return f"agostic {parent_sym}-H → σ-complex dative donation to M"
                # Dihydrogen (single H of pair — find partner)
                for other_h in adj_lewis[rep]:
                    if atoms[other_h] == 'H':
                        return "η²-H₂ dihydrogen — σ(H-H) dative donor"
                return "H σ-complex — dative donor"

        if cbc_char == 'L':
            is_CO = (sym == 'C' and any(atoms[k] == 'O' and bo_ij(rep, k) >= 2
                                        for k in adj_lewis[rep]))
            if is_CO:
                return "CO — σ-donor via C lone pair"
            if sym in ('C', 'N') and n_sub <= 2:
                return f"{sym} carbene/NHC — lone-pair σ-donor"
            # Anionic oxyanion O with 1 bond (NO3-, CO3-, SO4-): dative via lone pair
            if sym == 'O' and n_sub < n_val and lp_cnt > 0:
                fc_rep = fc[rep] if fc else 0
                if fc_rep < 0:
                    return f"O⁻ (oxyanion) — lone-pair dative donor ({lp_cnt} LP)"
                return f"O=N/S (resonance) — lone-pair dative donor ({lp_cnt} LP)"
            if lp_cnt > 0:
                return f"{sym}: {n_sub}/{n_val} bonds → lone-pair donor ({lp_cnt} LP)"
            return f"{sym}: {n_sub}/{n_val} bonds → dative donor"

        if cbc_char == 'X':
            # Use n_sub_total (bonds to all non-TM atoms) for an accurate description
            # when the atom is bonded to inorganic heteroatoms (S, P, Si, etc.)
            n_sub_total = sum(bo_ij(rep, k) for k in adj_lewis[rep]
                              if not is_TM(atoms[k]))
            deficit_real = n_val - n_sub_total
            if deficit_real <= 0:
                # Atom has its normal bonds but still classified X (e.g. sigma-bonded
                # to M in Lewis, or anionic with M-L covalent bond)
                return f"{sym}: {n_sub_total}/{n_val} bonds + covalent σ bond to M"
            return f"{sym}: {n_sub_total}/{n_val} bonds → {deficit_real} short → covalent σ to M"
        return ""

    for metal_idx in sorted(results):
        sym_m = atoms[metal_idx]
        interaction_records = results[metal_idx]  # list of (atom_tuple, cbc_char)

        total_L = sum(1 for _, t in interaction_records if t == 'L')
        total_X = sum(1 for _, t in interaction_records if t == 'X')
        total_Z = sum(1 for _, t in interaction_records if t == 'Z')

        parts = []
        if total_L: parts.append(f"L{_subscript(total_L)}")
        if total_X: parts.append(f"X{_subscript(total_X)}")
        if total_Z: parts.append(f"Z{_subscript(total_Z)}")
        designation = "".join(parts) if parts else "—"

        print(f"\n  {sym_m}{metal_idx+1}   [{designation}]")
        print(f"  {'Neighbour':<16} {'Avg.Dist':>9}  {'CBC':<5}  Notes")
        print(f"  {'-'*65}")

        # Collapse repeated (same atom_tuple, same cbc_char) rows into one with ×n
        from itertools import groupby
        seen_multi = {}   # (atom_tuple, cbc_char) -> count
        collapsed = []
        for atom_tuple, cbc_char in interaction_records:
            key = (atom_tuple, cbc_char)
            if key in seen_multi:
                seen_multi[key] += 1
            else:
                seen_multi[key] = 1
                collapsed.append((atom_tuple, cbc_char))

        for atom_tuple, cbc_char in collapsed:
            key = (atom_tuple, cbc_char)
            count = seen_multi[key]
            lbl   = row_label(atom_tuple)
            if count > 1:
                lbl = f"{lbl}(×{count})"
            d_str = f"{avg_dist(metal_idx, atom_tuple):.2f} Å"
            note  = row_note(atom_tuple, cbc_char, metal_idx)
            if count > 1:
                # Amend note for multi-bond terminal atoms
                sym_rep = atoms[atom_tuple[0]]
                bond_names = {2: 'double (2×)', 3: 'triple (3×)'}
                bond_desc = bond_names.get(count, f'{count}×')
                note = f"{sym_rep}: {bond_desc} covalent M={sym_rep} bond (nitrido/oxo/imido)"
            print(f"  {lbl:<16} {d_str:>9}  {cbc_char:<5}  {note}")

        print(f"\n  Total:  {total_L}L + {total_X}X + {total_Z}Z  →  {designation}")


# ─── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("="*62)
        print(" NBO Best Lewis Structure + $CHOOSE")
        print(" Demo: methylamine CH3NH2  (coords from nbo.log)")
        print("="*62)
        atoms, coords, charge = DEMO_ATOMS, DEMO_COORDS, 0
    elif sys.argv[1].endswith('.xyz'):
        atoms, coords, charge = read_xyz(sys.argv[1])
        charge = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        print(f"Read {len(atoms)} atoms from {sys.argv[1]}  (charge={charge})")
    elif sys.argv[1].endswith('.47'):
        atoms, coords, charge = read_nbo47(sys.argv[1])
        charge = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        print(f"Read {len(atoms)} atoms from {sys.argv[1]}  (charge={charge})")
    else:
        print(__doc__); sys.exit(0)

    unknown = sorted(set(a for a in atoms if a not in VALENCE))
    if unknown:
        print(f"WARNING: Unknown elements {unknown} – treated as C-like (val=4)")
        for u in unknown:
            VALENCE[u]=4; CORE_E[u]=0; COV_R[u]=0.77; STD_CAP[u]=4; ENEG[u]=2.55

    bo, lp, fc, stats = find_best_lewis(atoms, coords, charge)
    print()
    print_summary_and_choose(atoms, bo, lp, stats)

    # ── Octet / valence satisfaction check ────────────────────────────────────
    adj_check = defaultdict(list)
    for (i, j) in bo:
        adj_check[i].append(j); adj_check[j].append(i)
    violations = check_octet_violations(atoms, bo, adj_check, lp=lp, fc=fc)
    if violations:
        print()
        print("=" * 72)
        print("  ⚠  Octet / Valence Violations in Best Lewis Structure")
        print("=" * 72)
        print("  Atoms where bond_e + 2×LP ≠ expected valence electrons.")
        print("  These indicate a suboptimal Kekulé pattern or unusual geometry.")
        print()
        print(f"  {'Atom':<10} {'Bond-e':>6}  {'LP':>4}  {'Total-e':>7}  {'Exp-e':>6}  {'FC':>4}  Notes")
        print(f"  {'-'*68}")
        for idx, label, be, lp_cnt, exp_val, fc_val in violations:
            total_e = be + 2 * lp_cnt
            adj_syms = ', '.join(f"{atoms[k]}{k+1}" for k in adj_check[idx])
            print(f"  {label:<10} {be:>6}  {lp_cnt:>4}  {total_e:>7}  {exp_val:>6}  {fc_val:>+4}  bonded to: {adj_syms}")
        print()
        print(f"  → {len(violations)} atom(s) violate expected valence in this Lewis structure.")
    else:
        print()
        print("  ✓  All atoms satisfy their expected valence (octet rule OK).")

    print_cbc_report(atoms, coords, bo, lp, fc, charge)
