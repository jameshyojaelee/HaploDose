"""Bounded counterexamples for deposited DNA identity and resolution operations."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from haplodose.provenance import start_run, finish_run, file_record


R_CHECK = r'''
args <- commandArgs(trailingOnly=TRUE)
root <- parse(args[1])
stopifnot(length(root)==1, as.character(root[[1]][[2]])=="baf_cnv_ginkgo")
fun <- root[[1]][[3]]
stopifnot(is.call(fun), identical(fun[[1]], as.name("function")))
expressions <- as.list(fun[[3]])[-1]
guards <- Filter(function(x) is.call(x) && identical(x[[1]],as.name("if")),expressions)
stopifnot(length(guards)==6)
env <- new.env(parent=baseenv())
env$rdr <- matrix(0,2,3,dimnames=list(c("f1","f2"),c("a","b","c")))
env$mbaf <- env$rdr; colnames(env$mbaf) <- c("a","c","b")
env$fixed_rdr <- env$rdr; env$binclusters <- env$rdr
messages <- character()
withCallingHandlers(eval(guards[[1]],env),message=function(m){messages <<- c(messages,conditionMessage(m));invokeRestart("muffleMessage")})
stopifnot(length(messages)==0, !identical(colnames(env$rdr),colnames(env$mbaf)))
env$mbaf <- env$rdr; colnames(env$mbaf) <- c("z","x","y")
withCallingHandlers(eval(guards[[1]],env),message=function(m){messages <<- c(messages,conditionMessage(m));invokeRestart("muffleMessage")})
stopifnot(length(messages)==1)
helpers <- Filter(function(x) is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name("seg2cluster")),expressions)
stopifnot(length(helpers)==1)
helper <- eval(helpers[[1]][[3]],env)
labels <- helper(matrix(rep(1,12),ncol=1,dimnames=list(paste0("f",1:12),"cell")))
stopifnot(identical(as.numeric(labels),rep(1,12)))
# Independent intended chromosome-contained runs for identical values on two
# chromosomes. The helper has no chromosome argument; this is an input limit.
resets <- c(TRUE,rep(FALSE,5),TRUE,rep(FALSE,5))
stopifnot(identical(as.numeric(cumsum(resets)),c(rep(1,6),rep(2,6))))
cat("partial_permutation_missed\t1\nfull_mismatch_message_only\t1\nchromosome_boundary_unrepresented\t1\n")
cat(R.version.string,"\n")
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    parent = args.task/'dna_author_source_inventory_v1'
    inputs = dict(script=Path(__file__), design=args.task/'DNA_SOURCE_COUNTEREXAMPLE_DESIGN.md',
        source_manifest=parent/'run_manifest.json')
    run = start_run(args.out, inputs, {}, dict(scope='Synthetic source counterexamples; no pilot result inference'),
                    mode='scone_dna_source_counterexamples')
    try:
        prior = json.loads((parent/'run_manifest.json').read_text()); assert prior['status'] == 'complete'
        for record in prior['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        binary = shutil.which('Rscript'); assert binary is not None
        script = args.out/'bounded_reference.R'; script.write_text(R_CHECK)
        result = subprocess.run([binary, '--vanilla', str(script), str(parent/'source/DNA_related/Calculate_integerCNV.R')],
            capture_output=True, text=True, timeout=120)
        (args.out/'R.stdout').write_text(result.stdout); (args.out/'R.stderr').write_text(result.stderr)
        result.check_returncode()
        assert all(f'{name}\t1' in result.stdout for name in ('partial_permutation_missed',
            'full_mismatch_message_only', 'chromosome_boundary_unrepresented'))
        # Independently calculate exact inclusive overlaps. A destination bin
        # straddles two adjacent coarse bins with discordant allele evidence.
        target = (4, 7); source = [(1, 5, .1), (6, 10, .4)]
        overlap = [max(0, min(target[1], end)-max(target[0], start)+1) for start,end,value in source]
        assert overlap == [2, 2] and sum(overlap) == 4
        first = [value for start,end,value in source if start <= target[1] and end >= target[0]][0]
        reversed_first = [value for start,end,value in reversed(source) if start <= target[1] and end >= target[0]][0]
        assert first == .1 and reversed_first == .4
        weighted_sum = sum(n*value for n,(_,_,value) in zip(overlap,source))
        assert weighted_sum == 1. and weighted_sum/sum(overlap) == .25
        # Discarded evidence is undetectable after first-overlap projection.
        changed = [(1,5,.1),(6,10,.49)]
        assert [v for a,b,v in changed if a <= target[1] and b >= target[0]][0] == first
        report = dict(Rscript=file_record(Path(binary).resolve()), R_output=result.stdout,
            exact_identity_guard_partial_permutation_missed=True, full_mismatch_stops=False,
            chromosome_boundary_not_represented=True,
            first_overlap_reference=dict(overlap_bases=overlap, first=first, reversed_first=reversed_first,
                length_weighted_descriptive_average=.25, discarded_change_invisible=True,
                executed_author_transfer_function=False),
            original_author_function_executed='Only parsed first ID guard and isolated seg2cluster helper on synthetic inputs',
            published_pilot_execution_verified=False, assay_values_parsed=False, DNA_CNA_truth_constructed=False,
            production_default_changed=False, biological_null_specimens_verified=0)
        (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
