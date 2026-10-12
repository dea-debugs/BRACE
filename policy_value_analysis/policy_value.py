"""Policy value and overlap diagnostics for BRACE (METABRIC test set, n=594).
Estimand: restricted mean survival time at horizon H under a treatment policy.
Estimators: IPW with IPCW (Hajek), and doubly robust (AIPW) with outcome models.
Nuisance models fitted on the 1,386 training patients (primary) and cross fitted on the test set (sensitivity).
"""
import pandas as pd, numpy as np, json, sys, warnings
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import KFold
warnings.filterwarnings('ignore')
rng = np.random.default_rng(2026)

U = '/mnt/user-data/uploads/4_Deanna/BRACE/Datasets/analysis_extracts/'
S = '/tmp/claude-0/-home-claude/b1f82c5d-573c-5bd8-b091-63beff174c30/scratchpad/pv/'
OUT = sys.argv[1] if len(sys.argv) > 1 else S
H = 150.0
PAM = ["ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C"]  # 49 used historically (UBE2T excluded)

# ---------------- data ----------------
test = pd.read_csv(U + 'metabric_test594_clin_pam50.csv')
allp = pd.read_csv(U + 'metabric_all1980_clin_pam50.csv')
pat = pd.read_csv(S + 'patient.txt', sep='\t', comment='#')
sam = pd.read_csv(S + 'sample.txt', sep='\t', comment='#')
clin = pat[['PATIENT_ID','HORMONE_THERAPY','RADIO_THERAPY','CLAUDIN_SUBTYPE']].merge(
    sam[['PATIENT_ID','GRADE','TUMOR_SIZE','TUMOR_STAGE','ER_STATUS','PR_STATUS','HER2_STATUS']], on='PATIENT_ID')

def prep(df, is_test):
    d = df.copy()
    if is_test:
        d = d.rename(columns={'time.x':'time','event.x':'event','treatment.x':'T'})
    else:
        d = d.rename(columns={'treatment':'T'})
    d = d.merge(clin, left_on='patientID', right_on='PATIENT_ID', how='left')
    d['HORMONE'] = (d['HORMONE_THERAPY'] == 'YES').astype(float)
    d['RADIO'] = (d['RADIO_THERAPY'] == 'YES').astype(float)
    d['HER2'] = (d['HER2_STATUS'] == 'Positive').astype(float)
    d['PR'] = (d['PR_STATUS'] == 'Positive').astype(float)
    return d

te = prep(test, True)
al = prep(allp, False)
tr = al[~al['patientID'].isin(te['patientID'])].copy()
assert len(tr) == 1386, len(tr)
te['R'] = te['recommendation'].astype(int)

CLIN = ['AGE_AT_DIAGNOSIS','LYMPH_NODES_EXAMINED_POSITIVE','NPI','ER_IHC','INFERRED_MENOPAUSAL_STATE','GRADE','TUMOR_SIZE','TUMOR_STAGE','HER2','PR']
COTX = ['HORMONE','RADIO']

def design(d, cols, ref=None):
    """Median imputation (from ref) plus missing indicators for columns with any missingness in ref."""
    ref = d if ref is None else ref
    X = pd.DataFrame(index=d.index)
    for c in cols:
        v = pd.to_numeric(d[c], errors='coerce')
        rv = pd.to_numeric(ref[c], errors='coerce')
        X[c] = v.fillna(rv.median())
        if rv.isna().any():
            X[c + '_miss'] = v.isna().astype(float)
    return X

SETS = {'clinical': CLIN, 'pam50': PAM, 'clinical+pam50': CLIN + PAM, 'clinical+pam50+cotreatment': CLIN + PAM + COTX}

# ---------------- survival helpers ----------------
def km(time, event, grid=None):
    """Kaplan Meier S(t) evaluated on a grid (right continuous)."""
    t = np.asarray(time, float); e = np.asarray(event, int)
    o = np.argsort(t); t, e = t[o], e[o]
    ut = np.unique(t[e == 1])
    S = np.ones(len(ut)); s = 1.0
    for k, u in enumerate(ut):
        n = (t >= u).sum(); d = ((t == u) & (e == 1)).sum()
        s *= 1 - d / n; S[k] = s
    def Sfun(x):
        x = np.atleast_1d(np.asarray(x, float))
        idx = np.searchsorted(ut, x, side='right') - 1
        return np.where(idx < 0, 1.0, S[np.clip(idx, 0, max(len(S)-1,0))]) if len(ut) else np.ones_like(x)
    return Sfun, ut, S

def rmst(time, event, h, weights=None):
    """RMST = integral_0^h S(t) dt, with optional case weights (weighted KM)."""
    t = np.asarray(time, float); e = np.asarray(event, int)
    w = np.ones(len(t)) if weights is None else np.asarray(weights, float)
    o = np.argsort(t); t, e, w = t[o], e[o], w[o]
    ut = np.unique(t[(e == 1) & (t <= h)])
    S = 1.0; area = 0.0; prev = 0.0
    for u in ut:
        area += S * (u - prev); prev = u
        n = w[t >= u].sum(); d = w[(t == u) & (e == 1)].sum()
        S *= 1 - d / n
    area += S * (h - prev)
    return area

def ipcw(d, h, by_arm=True):
    """Weights for the restricted outcome: Delta_h / G(min(U,h)-), G = KM of censoring, by treatment arm."""
    U_ = np.minimum(d['time'].values, h)
    obs = ((d['event'].values == 1) & (d['time'].values <= h)) | (d['time'].values >= h)
    G = np.ones(len(d))
    groups = [d['T'].values == 1, d['T'].values == 0] if by_arm else [np.ones(len(d), bool)]
    for g in groups:
        Sfun, _, _ = km(d['time'].values[g], 1 - d['event'].values[g])
        # G evaluated just before min(U,h): use U- by subtracting a tiny epsilon
        G[g] = Sfun(U_[g] - 1e-9)
    G = np.clip(G, 0.05, 1.0)
    return U_, obs.astype(float), G

# ---------------- nuisance models ----------------
PS_C = 0.2
def fit_propensity(Xtr, Ttr):
    m = make_pipeline(StandardScaler(), LogisticRegression(C=PS_C, max_iter=2000))
    m.fit(Xtr, Ttr); return m

def fit_outcome(Xtr, Ttr, Y, w):
    """m_t(x) = E[min(S,h) | X, T=t], ridge on IPCW weighted observed restricted times."""
    models = {}
    for t in (0, 1):
        idx = (Ttr == t) & (w > 0)
        m = make_pipeline(StandardScaler(), Ridge(alpha=10.0))
        m.fit(Xtr[idx], Y[idx], ridge__sample_weight=w[idx]); models[t] = m
    return models

def policy_values(te_, e_hat, m0, m1, pi_dict, h, trim=(0.02, 0.98)):
    U_, obs, G = ipcw(te_, h)
    T = te_['T'].values; Yh = U_
    e_hat = np.clip(e_hat, *trim)
    cw = obs / G  # censoring weight
    res = {}
    for name, pi in pi_dict.items():
        pi = np.asarray(pi, int)
        ep = np.where(pi == 1, e_hat, 1 - e_hat)
        ind = (T == pi).astype(float)
        w = ind * cw / ep
        ipw = np.sum(w * Yh) / np.sum(w)                      # Hajek
        mpi = np.where(pi == 1, m1, m0)
        dr = np.mean(mpi + ind * cw / ep * (Yh - mpi))
        res[name] = {'IPW': ipw, 'DR': dr, 'ESS': w.sum()**2 / np.sum(w**2), 'n_matching': int(ind.sum())}
    # observed policy: plain IPCW mean (equals KM RMST)
    res['observed'] = {'IPW': np.sum(cw * Yh) / np.sum(cw), 'DR': np.sum(cw * Yh) / np.sum(cw), 'ESS': cw.sum()**2 / np.sum(cw**2), 'n_matching': len(T)}
    return res

def run(adj, nuisance='train', seed=None, te_=None, tr_=None, overlap=None):
    te_ = te if te_ is None else te_; tr_ = tr if tr_ is None else tr_
    cols = SETS[adj]
    Xte = design(te_, cols, ref=tr_).values
    pis = {'BRACE': te_['R'].values, 'treat_none': np.zeros(len(te_), int), 'treat_all': np.ones(len(te_), int)}
    if nuisance == 'train':
        Xtr = design(tr_, cols).values
        Utr, obstr, Gtr = ipcw(tr_, H)
        ps = fit_propensity(Xtr, tr_['T'].values)
        om = fit_outcome(Xtr, tr_['T'].values, Utr, obstr / Gtr)
        e_hat = ps.predict_proba(Xte)[:, 1]; m0 = om[0].predict(Xte); m1 = om[1].predict(Xte)
        if overlap is not None:  # restrict evaluation to the overlap population
            keep = (e_hat >= overlap[0]) & (e_hat <= overlap[1])
            te_ = te_[keep].reset_index(drop=True); e_hat, m0, m1 = e_hat[keep], m0[keep], m1[keep]
            pis = {k: v[keep] for k, v in pis.items()}
            pv = policy_values(te_, e_hat, m0, m1, pis, H); pv['n_overlap'] = int(keep.sum())
            pv['reclassified_in_overlap'] = int(((te_['T'] != te_['R'])).sum())
            return pv, e_hat
    else:  # cross fit on test
        e_hat = np.zeros(len(te_)); m0 = np.zeros(len(te_)); m1 = np.zeros(len(te_))
        Ute, obste, Gte = ipcw(te_, H)
        for trn, val in KFold(5, shuffle=True, random_state=seed or 1).split(Xte):
            ps = fit_propensity(Xte[trn], te_['T'].values[trn])
            om = fit_outcome(Xte[trn], te_['T'].values[trn], Ute[trn], obste[trn] / Gte[trn])
            e_hat[val] = ps.predict_proba(Xte[val])[:, 1]; m0[val] = om[0].predict(Xte[val]); m1[val] = om[1].predict(Xte[val])
    return policy_values(te_, e_hat, m0, m1, pis, H), e_hat

# ---------------- 0. validation against the paper ----------------
val = {}
for arm in (0, 1):
    d = te[te['T'] == arm]
    f = d[d['followed'] == 1]; nf = d[d['followed'] == 0]
    val[f'arm{arm}'] = rmst(f['time'], f['event'], H) - rmst(nf['time'], nf['event'], H)
f = te[te['followed'] == 1]; nf = te[te['followed'] == 0]
val['pooled'] = rmst(f['time'], f['event'], H) - rmst(nf['time'], nf['event'], H)
print('Validation (paper: no chemo arm 19.19, chemo arm 6.24, pooled 18.94):', {k: round(v, 2) for k, v in val.items()})

# ---------------- 1. diagnostics ----------------
diag = {}
Utr, obstr, Gtr = ipcw(tr, H)
for adj in SETS:
    cols = SETS[adj]
    Xtr = design(tr, cols).values; Xte = design(te, cols, ref=tr).values
    ps = fit_propensity(Xtr, tr['T'].values); e = ps.predict_proba(Xte)[:, 1]
    T = te['T'].values
    q = lambda v: np.percentile(v, [1, 5, 25, 50, 75, 95, 99]).round(3).tolist()
    w = np.where(T == 1, 1 / e, 1 / (1 - e))
    # standardised mean differences on clinical covariates, before and after IPW
    Xc = design(te, CLIN, ref=tr)
    smd = {}
    for c in Xc.columns:
        x = Xc[c].values
        m1_, m0_ = x[T == 1].mean(), x[T == 0].mean(); s = np.sqrt((x[T == 1].var() + x[T == 0].var()) / 2) + 1e-9
        wm1 = np.average(x[T == 1], weights=w[T == 1]); wm0 = np.average(x[T == 0], weights=w[T == 0])
        smd[c] = {'before': round((m1_ - m0_) / s, 3), 'after_IPW': round((wm1 - wm0) / s, 3)}
    diag[adj] = {
        'propensity_treated_pct': q(e[T == 1]), 'propensity_untreated_pct': q(e[T == 0]),
        'share_outside_0.05_0.95': round(float(((e < 0.05) | (e > 0.95)).mean()), 3),
        'ESS_treated': round(float(w[T == 1].sum()**2 / (w[T == 1]**2).sum()), 1), 'n_treated': int((T == 1).sum()),
        'ESS_untreated': round(float(w[T == 0].sum()**2 / (w[T == 0]**2).sum()), 1), 'n_untreated': int((T == 0).sum()),
        'propensity_among_BRACE_rec_chemo_pct': q(e[te['R'].values == 1]),
        'smd_clinical': smd,
    }

# missingness and follow up
def rev_km_median(time, event):
    Sfun, ut, S = km(time, 1 - np.asarray(event))
    idx = np.where(S <= 0.5)[0]
    return float(ut[idx[0]]) if len(idx) else float('nan')
miss = {}
for name, d in (('test594', te), ('train1386', tr), ('all1980', al)):
    miss[name] = {
        'n': len(d), 'deaths': int(d['event'].sum()), 'median_follow_up_reverse_KM': round(rev_km_median(d['time'], d['event']), 1),
        'share_followed_150m_or_event': round(float(((d['time'] >= H) | (d['event'] == 1)).mean()), 3),
        'censored_before_150m': int(((d['time'] < H) & (d['event'] == 0)).sum()),
        'missing': {c: int(pd.to_numeric(d[c], errors='coerce').isna().sum()) for c in CLIN + COTX},
    }

# ---------------- 2. policy values ----------------
results = {}
for adj in SETS:
    for nuis in ('train', 'crossfit'):
        pv, _ = run(adj, nuis)
        results[f'{adj}|{nuis}'] = pv
    pv, _ = run(adj, 'train', overlap=(0.1, 0.9)); results[f'{adj}|train|overlap_0.1_0.9'] = pv

# ---------------- 3. bootstrap CIs (primary: clinical+pam50, nuisance on train) ----------------
B = int(sys.argv[2]) if len(sys.argv) > 2 else 300
boot = {k: [] for k in ['BRACE_IPW', 'BRACE_DR', 'obs', 'none_IPW', 'none_DR', 'all_IPW', 'all_DR', 'ov_BRACE_DR', 'ov_none_DR', 'ov_all_DR', 'ov_obs', 'ov_BRACE_IPW', 'ov_none_IPW']}
for b in range(B):
    teb = te.sample(len(te), replace=True, random_state=b).reset_index(drop=True)
    trb = tr.sample(len(tr), replace=True, random_state=10_000 + b).reset_index(drop=True)
    try:
        pv, _ = run('clinical+pam50', 'train', te_=teb, tr_=trb)
    except Exception as ex:
        continue
    boot['BRACE_IPW'].append(pv['BRACE']['IPW']); boot['BRACE_DR'].append(pv['BRACE']['DR']); boot['obs'].append(pv['observed']['IPW'])
    boot['none_IPW'].append(pv['treat_none']['IPW']); boot['none_DR'].append(pv['treat_none']['DR'])
    boot['all_IPW'].append(pv['treat_all']['IPW']); boot['all_DR'].append(pv['treat_all']['DR'])
    pvo, _ = run('clinical+pam50', 'train', te_=teb, tr_=trb, overlap=(0.1, 0.9))
    boot['ov_BRACE_DR'].append(pvo['BRACE']['DR']); boot['ov_none_DR'].append(pvo['treat_none']['DR']); boot['ov_all_DR'].append(pvo['treat_all']['DR']); boot['ov_obs'].append(pvo['observed']['IPW'])
    boot['ov_BRACE_IPW'].append(pvo['BRACE']['IPW']); boot['ov_none_IPW'].append(pvo['treat_none']['IPW'])
bt = {k: np.array(v) for k, v in boot.items()}
def ci(a): return [round(float(np.percentile(a, 2.5)), 2), round(float(np.percentile(a, 97.5)), 2)]
ci_out = {
    'B': len(bt['obs']),
    'V(BRACE) IPW': ci(bt['BRACE_IPW']), 'V(BRACE) DR': ci(bt['BRACE_DR']), 'V(observed)': ci(bt['obs']),
    'V(BRACE) - V(observed) IPW': ci(bt['BRACE_IPW'] - bt['obs']), 'V(BRACE) - V(observed) DR': ci(bt['BRACE_DR'] - bt['obs']),
    'V(BRACE) - V(treat none) IPW': ci(bt['BRACE_IPW'] - bt['none_IPW']), 'V(BRACE) - V(treat none) DR': ci(bt['BRACE_DR'] - bt['none_DR']),
    'V(treat all) - V(treat none) DR': ci(bt['all_DR'] - bt['none_DR']),
    'overlap: V(BRACE) - V(observed) DR': ci(bt['ov_BRACE_DR'] - bt['ov_obs']),
    'overlap: V(BRACE) - V(treat none) DR': ci(bt['ov_BRACE_DR'] - bt['ov_none_DR']),
    'overlap: V(BRACE) - V(treat none) IPW': ci(bt['ov_BRACE_IPW'] - bt['ov_none_IPW']),
    'overlap: V(treat all) - V(treat none) DR': ci(bt['ov_all_DR'] - bt['ov_none_DR']),
}

json.dump({'validation': val, 'diagnostics': diag, 'followup_missingness': miss, 'policy_values': results, 'bootstrap_CI': ci_out},
          open(OUT + 'policy_value_results.json', 'w'), indent=1, default=float)
print(json.dumps({'policy_values': results, 'bootstrap_CI': ci_out}, indent=1, default=float))
