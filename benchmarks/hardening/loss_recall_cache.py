"""Research-only exact lookup of strict-threshold loss recall, not confidence."""
import numpy as np


class RecallCurve:
    def __init__(self, scores, loss, common):
        scores = np.asarray(scores, float)
        loss, common = np.asarray(loss, bool), np.asarray(common, bool)
        if scores.ndim != 2 or loss.shape != scores.shape or common.shape != scores.shape:
            raise ValueError('Aligned window by patient arrays required')
        mask = loss & common
        self.denominator = mask.sum(axis=0)
        # NaN scores remain in the denominator, as in the original scorer;
        # their comparisons are always false. +/-inf retain comparison semantics.
        values = [np.sort(scores[:, j][mask[:, j] & ~np.isnan(scores[:, j])])
                  for j in range(scores.shape[1])]
        self.support = np.unique(np.concatenate(values)) if values else np.empty(0)
        self.counts = np.empty((len(self.support) + 1, scores.shape[1]), np.int64)
        for j, v in enumerate(values):
            self.counts[:-1, j] = np.searchsorted(v, self.support, side='left')
            self.counts[-1, j] = len(v)
        for array in (self.denominator, self.support, self.counts):
            array.setflags(write=False)

    def at(self, threshold):
        result = np.full(len(self.denominator), np.nan)
        if np.isfinite(threshold):
            row = int(np.searchsorted(self.support, threshold, side='left'))
            np.divide(self.counts[row], self.denominator, out=result, where=self.denominator > 0)
        return result
