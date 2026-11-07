"""Compute device for the JAX kernels: CPU (default), GPU (opt in), or auto.

JAX binds to a platform the first time it computes anything, so the device must be chosen before that:
`haplodose run/bulk --device` and `run(..., device=)` / `run_bulk(..., device=)` call `select` first. Importing
haplodose does not start JAX.

- 'cpu': the default. Reproducible across machines with the same CPU generation; never touches CUDA.
- 'gpu': CUDA through JAX (install the extra: `pip install "haplodose[cuda12]"`). Fails with a clear error when no
  usable GPU is found instead of falling back silently. The HMM switches to the associative scan on GPU
  (`kernels.hmm.default_method`), and float reduction order differs, so borderline calls can move slightly.
- 'auto': GPU when JAX finds one, else CPU.
Only the JAX kernels (HMM, likelihoods, cell E-step, programme scores) run on the device; Leiden (igraph), PCA
(scikit-learn) and the clone search bookkeeping stay on the CPU.
"""
from __future__ import annotations

DEVICES = ('cpu', 'gpu', 'auto')
_PLATFORM = {'cpu': 'cpu', 'gpu': 'cuda', 'auto': None}


def _started() -> bool:
    from jax._src import xla_bridge
    return xla_bridge.backends_are_initialized()


def select(device: str = 'cpu') -> str:
    """Bind JAX to `device` and return the backend in use ('cpu' or 'gpu')."""
    if device not in DEVICES:
        raise ValueError(f'device must be one of {DEVICES}, got {device!r}')
    import jax
    if _started():
        backend = jax.default_backend()
        if device == 'auto' or backend == device:
            return backend
        raise RuntimeError(f'JAX is already running on {backend}; choose the device before any JAX computation '
                           f'(call haplodose.device.select({device!r}) first, or start a new process)')
    jax.config.update('jax_platforms', _PLATFORM[device])
    try:
        backend = jax.default_backend()
    except Exception as err:   # no CUDA: RuntimeError, or AssertionError from the cuda12 plugin (jax 0.10.2)
        if device == 'gpu':
            raise RuntimeError('device gpu: JAX found no usable CUDA GPU. Install the CUDA extra '
                               '(pip install "haplodose[cuda12]"), check that a GPU is visible (nvidia-smi), '
                               'or use --device cpu.') from err
        raise
    if device == 'gpu' and backend != 'gpu':
        raise RuntimeError(f'device gpu: JAX started on {backend}, not a GPU; use --device cpu or auto')
    return backend
