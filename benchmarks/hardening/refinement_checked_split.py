"""Research split stability check using the unchanged EM and merge criteria.

This adds a search constraint, not calibration. Preview computations never
replace accepted inference states. Nonspatial expression-only development.
"""
from contextlib import contextmanager
from copy import deepcopy
import numpy as np
from haplodose.clones import SpecimenModel


def snapshot(model):
    # Rows is a read-only data accessor: preserve its identity for the frozen
    # gene-fold registry. Its fold draws are deterministic functions of row IDs.
    return deepcopy(model.__dict__, {id(model):model, id(model.rows):model.rows,
                                    id(model.params):model.params})


def restore(model, state):
    model.__dict__.clear()
    model.__dict__.update(state)


@contextmanager
def checked_splits(records, *, apply):
    original = SpecimenModel.try_splits

    def checked(self, R, doses, pi):
        if self.alleles is not None or self.graph is not None:
            raise ValueError('Refinement-check research supports expression-only nonspatial fits')
        before = snapshot(self)
        old_doses, old_pi = [np.asarray(d).copy() for d in doses], pi.copy()
        result = original(self, R, doses, pi)
        if not result[2]:
            return result
        proposed = snapshot(self)
        log_start = len(self.log)
        trial_doses, trial_pi = result[:2]
        rounds = []
        try:
            # A merge strictly reduces the finite component count. Apply the
            # existing complete refit/merge check without proposing more groups.
            for _ in range(len(trial_doses)):
                fitted = self.em(trial_doses, trial_pi)
                if len(fitted[0]) > len(trial_doses):
                    raise RuntimeError('Preview EM added a group')
                trial_doses, trial_pi, merged = self.merge_unsupported(fitted[2], fitted[0], fitted[1])
                rounds.append(dict(fitted_groups=len(fitted[0])-1,
                                   remaining_groups=len(trial_doses)-1, merged=bool(merged)))
                if not merged:
                    break
                if len(trial_doses) >= len(fitted[0]):
                    raise RuntimeError('Preview merge did not reduce group count')
            else:
                raise RuntimeError('Split preview did not terminate')
            passed = len(trial_doses) > len(old_doses)
            entry = dict(before_groups=len(old_doses)-1, proposed_groups=len(result[0])-1,
                         preview_groups=len(trial_doses)-1, passed=passed, apply=apply,
                         rounds=rounds, preview_log=deepcopy(self.log[log_start:]),
                         criteria_changed=False, calibration_claim=False)
            records.append(entry)
        except BaseException:
            restore(self, proposed)
            raise
        # On acceptance or sham, the original post-proposal state is restored.
        # The outer controller executes its usual EM, without extra fit rounds.
        if not apply or passed:
            restore(self, proposed)
            return result
        history = proposed['log']
        restore(self, before)
        self.log = history
        self.log.append(dict(action='reject_split_after_refinement',
                             before_groups=entry['before_groups'],
                             proposed_groups=entry['proposed_groups'],
                             preview_groups=entry['preview_groups'],
                             reason='no_supported_additional_group_after_unchanged_refit_and_merge'))
        # No accepted addition: stop this search rather than repeatedly proposing
        # the same failed addition. Original doses, pi, R/LL and nuisances agree.
        return old_doses, old_pi, False

    SpecimenModel.try_splits = checked
    try:
        yield
    finally:
        SpecimenModel.try_splits = original
