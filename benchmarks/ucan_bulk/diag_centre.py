"""Does Haplodose bulk's per-sample centre (fit_pooled: weighted median of residuals after the HMM dose) cost
across-patient arm AUC? Recomputes arm medians from a run's log_ratio.npy with the centre added back and the
plain median over bins subtracted instead. Usage: diag_centre.py BULK_RUN OUT_ARMS_TSV
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

run, out = Path(sys.argv[1]), sys.argv[2]
r = np.load(run / 'log_ratio.npy').astype(np.float64)
s = pd.read_csv(run / 'samples.tsv', sep='\t')
arm = pd.read_csv(run / 'bins.tsv', sep='\t').arm.to_numpy()
r = r + s.centre.to_numpy()[:, None]
r -= np.median(r, axis=1, keepdims=True)
A = pd.DataFrame(r.T, columns=s['sample']).groupby(arm).median().T
A.index.name = 'sample'
A.reset_index().melt(id_vars='sample', var_name='arm', value_name='value').to_csv(out, sep='\t', index=False)
