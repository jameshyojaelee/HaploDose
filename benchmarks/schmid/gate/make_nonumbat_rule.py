"""Append a no-Numbat copy of one of their evaluation rules to a WORKING COPY of workflow.sm (never the pinned clone).

The copy is the rule text verbatim with three literal edits: the rule name gets '_nonumbat', the numbat_obj input
line is removed (their evaluate_results.R skips Numbat when snakemake@input$numbat_obj is NULL, lines 193-217), and
'evaluation/' and the log path are redirected to 'evaluation_nonumbat/'. Numbat bins then do not enter the combined
bin intersection, so these values are provisional gate rows (STATUS claim 5, 2026-11-04).
Usage: make_nonumbat_rule.py WORKFLOW RULE_NAME
"""
import re
import sys


def main(path, name):
    s = open(path).read()
    m = re.search(rf'^rule {re.escape(name)}:\n.*?(?=^rule |\Z)', s, flags=re.M | re.S)
    assert m, name
    r = m.group(0)
    n_numbat = len(re.findall(r'^\s*numbat_obj = .*\n', r, flags=re.M))
    assert n_numbat == 1, n_numbat
    r = re.sub(r'^\s*numbat_obj = .*\n', '', r, flags=re.M)
    assert r.count('/evaluation/') >= 5
    r = r.replace('/evaluation/', '/evaluation_nonumbat/')
    r = r.replace(f'rule {name}:', f'rule {name}_nonumbat:', 1)
    r = re.sub(r'(\t\t")logs/', r'\1logs/nonumbat_', r)
    assert 'numbat' not in r.replace('_nonumbat', '').replace('nonumbat_', '')
    with open(path, 'a') as f:
        f.write('\n\n# Added for the claim 5 provisional gate (no Numbat); see make_nonumbat_rule.py\n' + r.rstrip() + '\n')
    print(r)


if __name__ == '__main__':
    main(*sys.argv[1:3])
