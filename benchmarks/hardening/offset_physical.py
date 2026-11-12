"""Constructor-compatible depth checks and strengthened tests of the built candidate."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    script = Path(__file__).with_name('offset_expectation.py')
    tests = Path(__file__).parents[2]/'tests/test_offset_conditional_candidate.py'
    rec = start_run(a.out, dict(script=Path(__file__), counter=script, tests=tests,
        design=a.task/'PHYSICAL_DEPTH_DESIGN.md', validation=a.task/'validation_23089312/pytest.txt'), {},
        dict(scope='Failed release candidate diagnostics, not relaxed validation'), mode='offset_physical')
    try:
        for mode, extra in [('before', ['--legacy-method']), ('after', ['--conditional-candidate', '--expect-fixed'])]:
            with (a.out/f'{mode}.txt').open('w') as handle:
                subprocess.run([sys.executable, '-I', str(script), '--out', str(a.out/mode),
                    '--task', str(a.task), '--physical-depth', *extra], check=True, stdout=handle, stderr=subprocess.STDOUT)
        with (a.out/'focused_pytest.txt').open('w') as handle:
            subprocess.run([sys.executable, '-I', '-m', 'pytest', '--basetemp', str(a.out/'pytest_tmp'),
                '--import-mode=importlib', '-p', 'no:cacheprovider', '-q', str(tests)],
                check=True, stdout=handle, stderr=subprocess.STDOUT)
        (a.out/'report.json').write_text(json.dumps(dict(physical_cases_per_arm=24,
            focused_tests_passed=True, release_validation_passed=False, original_recovery_assertion_preserved=True,
            biological_specimens=0, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU constructor-compatible exposures, archived and installed methods, strengthened pytest')


if __name__ == '__main__':
    main()
