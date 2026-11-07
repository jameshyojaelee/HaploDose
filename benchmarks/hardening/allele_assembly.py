"""Bounded assembly forensics for the evidence-selected allele trace; run on compute."""
from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

from haplodose.io import new_output_dir, write_json
from haplodose.provenance import file_record, start_run, finish_run


def text_head(path, limit=100):
    opener = gzip.open if str(path).endswith('.gz') else open
    lines = []
    with opener(path, 'rt') as handle:
        for _, line in zip(range(limit), handle):
            lines.append(line.rstrip())
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root', required=True, type=Path)
    ap.add_argument('--trace-root', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--bcftools', type=Path, help='Existing bcftools executable; otherwise conda then PATH')
    ap.add_argument('--template-only', type=Path, metavar='PREVIOUS_ROOT',
                    help='Capture installed Numbat template only, linking previous assembly evidence')
    args = ap.parse_args()
    root, trace = args.data_root.resolve(), args.trace_root.resolve()
    if args.out.exists():
        raise FileExistsError(f'Refusing existing output: {args.out}')
    out = new_output_dir(args.out)
    summary = json.loads((trace / 'summary.json').read_text())
    source = Path(__file__).resolve().parents[2]
    donor = root / 'results/kinker_numbat_20261102/pool16/alleles' / summary['selected_line']
    panel = root / 'data/external/numbat_hg19_panel'
    env = root.parent / 'CRC/data/interim/numbat_env'
    inputs = {'trace_manifest': trace / 'run_manifest.json', 'trace_summary': trace / 'summary.json',
              'runner': Path(__file__), 'phase_log': donor / 'phasing.log'}
    if args.template_only is not None:
        inputs['previous_assembly_manifest'] = args.template_only.resolve() / 'run_manifest.json'
        inputs['previous_assembly_commands'] = args.template_only.resolve() / 'commands.json'
        inputs['previous_assembly_summary'] = args.template_only.resolve() / 'summary.json'
    for name in ['phase_pool.R', 'pileup_channel.sbatch', 'build_hg19_panel.sbatch', 'numbat_pool.sbatch']:
        inputs[name] = source / 'benchmarks/kinker/alleles' / name
    for name in ['DESCRIPTION', 'data/Rdata.rdb', 'data/Rdata.rdx', 'data/Rdata.rds']:
        path = env / 'lib/R/library/numbat' / name
        if path.is_file():
            inputs['installed_numbat_' + name] = path
    rec = start_run(out, inputs, {}, {'purpose': 'assembly_disagreement_forensics',
        'limits': 'Headers and two target-defined indexed loci; no full BAM hash or reprocessing'}, mode='allele_assembly')
    commands, gaps = [], []

    def command(label, argv, timeout=45):
        entry = dict(label=label, argv=list(map(str, argv)))
        try:
            p = subprocess.run(entry['argv'], capture_output=True, text=True, timeout=timeout)
            entry.update(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)
            if p.returncode:
                gaps.append(label + ': command failed')
        except (OSError, subprocess.TimeoutExpired) as error:
            entry['error'] = str(error)
            gaps.append(label + ': unavailable or timed out')
        commands.append(entry)
        return entry.get('stdout', '')

    try:
        tools = {name: str(env / 'bin' / name) if (env / 'bin' / name).is_file() else shutil.which(name)
                 for name in ['samtools', 'bcftools', 'Rscript']}
        if args.bcftools is not None:
            tools['bcftools'] = str(args.bcftools.resolve())
        tool_records = {k: file_record(v) for k, v in tools.items() if v}
        for name in ['samtools', 'bcftools']:
            if tools[name]:
                command(name + '_version', [tools[name], '--version'])
        write_json(dict(LOADEDMODULES=os.environ.get('LOADEDMODULES'),
            MODULEPATH=os.environ.get('MODULEPATH'), tools=tool_records), out / 'tool_environment.json')
        template_expr = ('library(numbat); cat("VERSION",as.character(packageVersion("numbat")),"\\n"); '
            'print(getFromNamespace("make_vcf_chr","numbat")); '
            'holder <- new.env(); data("vcf_meta",package="numbat",envir=holder); '
            'stopifnot(exists("vcf_meta",envir=holder,inherits=FALSE)); '
            'cat("VCF_META_TEMPLATE\\n"); print(holder$vcf_meta); sessionInfo()')
        if args.template_only is not None:
            if not tools['Rscript']:
                raise RuntimeError('Rscript unavailable for template capture')
            command('installed_template_and_header_generator', [tools['Rscript'], '-e', template_expr], timeout=60)
            write_json(commands, out / 'commands.json')
            write_json(dict(previous_assembly_root=str(args.template_only.resolve()), gaps=gaps,
                purpose='Template-only follow-up; paired reference loci and BAM headers remain in linked previous execution'), out / 'summary.json')
            if gaps:
                raise RuntimeError('Template capture failed; see commands.json')
            finish_run(out, rec, backend='cpu')
            return
        phase_log = (donor / 'phasing.log').read_text()
        resolved = sorted(set(re.findall(r'(?:data/[^\s\\]+|/gpfs/[^\s\\]+)', phase_log)))
        write_json(dict(resolved_paths=[dict(recorded=p, resolved=str((root / p).resolve()),
            exists=(root / p).exists()) for p in resolved], excerpt=text_head(donor / 'phasing.log', 70)), out / 'phase_log_context.json')
        target = donor / 'phasing' / (summary['selected_line'] + '_chr1.vcf.gz')
        header, loci = [], []
        with gzip.open(target, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    header.append(line.rstrip())
                elif len(loci) < 2:
                    fields = line.split('\t')
                    loci.append(dict(chrom=fields[0], pos=int(fields[1]), ref=fields[3], alt=fields[4]))
                if len(loci) == 2:
                    break
        write_json(dict(header=header, loci=loci, selection='First two chr1 genotype-target rows, independent of fitted clone results'), out / 'target.json')
        phased = donor / 'phasing' / (summary['selected_line'] + '_chr1.phased.vcf.gz')
        write_json(dict(file=file_record(phased), header=[x for x in text_head(phased, 100) if x.startswith('#')]), out / 'phased_header.json')
        for name, path in [('panel', panel / 'panel/chr1.genotypes.bcf'),
                           ('sites', panel / 'snps_AF5e2_chr1to22_hg19.vcf.gz')]:
            if tools['bcftools']:
                command(name + '_header', [tools['bcftools'], 'view', '-h', path])
            else:
                gaps.append(name + ': bcftools unavailable')
        map_path = panel / 'genetic_map_hg19_withX.txt.gz'
        first, last, count = [], None, 0
        with gzip.open(map_path, 'rt') as f:
            for line in f:
                count += 1
                if count <= 4:
                    first.append(line.rstrip())
                last = line.rstrip()
        write_json(dict(file=file_record(map_path), first=first, last=last, lines=count,
            assembly_evidence='Executed Eagle path and panel construction recipe; map filename alone is not assembly validation'), out / 'genetic_map.json')
        bam_lengths, locus_checks = {}, []
        for index, (channel, assoc) in enumerate(sorted(summary['source_bam_associations'].items())):
            log = root / 'results/kinker_numbat_20261102/pileup' / channel / 'cellsnp_time.txt'
            write_json(dict(file=file_record(log), excerpt=text_head(log, 35)), out / (channel + '_pileup_log.json'))
            if not tools['samtools']:
                gaps.append(channel + ': samtools unavailable')
                continue
            bam = Path(assoc['bam'])
            h = command(channel + '_bam_header', [tools['samtools'], 'view', '-H', bam])
            bam_lengths[channel] = [line for line in h.splitlines() if line.startswith('@SQ') and re.search(r'\tSN:(?:chr)?1\t', line)]
            if index < len(loci):
                locus = loci[index]
                chrom = 'chr1' if any('\tSN:chr1\t' in x for x in bam_lengths[channel]) else '1'
                region = f'{chrom}:{locus["pos"]}-{locus["pos"]}'
                command(channel + '_single_locus_count', [tools['samtools'], 'view', '-c', bam, region])
                if tools['bcftools']:
                    observed = command(channel + '_panel_locus', [tools['bcftools'], 'query', '-r',
                        f'chr1:{locus["pos"]}-{locus["pos"]}', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', panel / 'panel/chr1.genotypes.bcf'])
                    expected = f'chr1\t{locus["pos"]}\t{locus["ref"]}\t{locus["alt"]}'
                    matches = expected in observed.splitlines()
                    locus_checks.append(dict(target=locus, panel_exact_allele_match=matches))
                    if not matches:
                        gaps.append(channel + ': target allele absent from indexed reference locus')
        if tools['Rscript']:
            # Inspect installed annotation and the exact generator of inherited VCF headers.
            expr = template_expr + '; print(utils::head(gtf_hg19)); print(utils::head(gtf_hg38))'
            command('installed_annotation_and_header_generator', [tools['Rscript'], '-e', expr], timeout=60)
        else:
            gaps.append('Rscript unavailable: installed annotation origin not checked')
        write_json(commands, out / 'commands.json')
        write_json(dict(tools=tool_records, bam_chr1_sequences=bam_lengths, locus_checks=locus_checks, gaps=gaps,
            interpretation='Compare actual BAM SQ lengths, panel/reference header provenance, executed map paths and Numbat header generator. Coordinate identity at two loci does not prove global assembly or donor identity.',
            known_chr1_lengths={'GRCh37': 249250621, 'GRCh38': 248956422},
            annotation_limit='R object previews and installed version recorded; source GTF accession and full gene/cM annotation accuracy remain unverified unless explicitly present in captured source.',
            bounded_loci=len(loci), bam_digest='not recomputed'), out / 'summary.json')
        finish_run(out, rec, backend='cpu')
    except Exception as error:
        write_json(commands, out / 'commands.json')
        finish_run(out, rec, backend='cpu', error=error)
        raise


if __name__ == '__main__':
    main()
