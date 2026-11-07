"""Download the MSK SPECTRUM DLP+ scDNA dataset (Synapse syn66366960, project syn66366718 "MSK SPECTRUM WGD").

Needs a Synapse login: anonymous users can list entities but not download files (checked 2026-11-02, the REST API
returns "Anonymous users have only READ access permission"). No access requirement (DUA) is attached to the project,
the dataset or a sampled file. One-time setup by the user:
    ! micromamba run -r /gpfs/commons/home/jameslee/micromamba -n synapse synapse config
(paste a Synapse personal access token with "view" and "download" scope; this writes ~/.synapseConfig).

The dataset has 239 files (benchmarks/heldout2/spectrum_synapse_entity_headers.json): per patient (41) signals_*.h5,
sbmclone_*_snv.h5; per patient (39) *_cell_info.csv.gz, *_cluster_info.csv.gz, *_branch_info.csv.gz,
*_clones_pruned.pickle; plus filtered_cell_table.csv.gz.
BLIND: download and checksum only; truth files are not opened beyond confirming they exist.
Writes <dest>/manifest.tsv (file, url, bytes, md5, md5_source, retrieved); Synapse verifies each download against
the file handle's contentMd5, and the md5 is recomputed locally here.

Usage: python -I fetch_spectrum_scdna_synapse.py <dest_dir>
"""
import datetime
import hashlib
import os
import sys

import synapseclient

DATASET = "syn66366960"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def main(dest):
    os.makedirs(dest, exist_ok=True)
    syn = synapseclient.Synapse(skip_checks=True)
    syn.login(silent=True)  # reads ~/.synapseConfig; fails here if no token
    dataset = syn.restGET(f"/entity/{DATASET}")
    items = [i["entityId"] for i in dataset["items"]]
    rows = []
    for k, ent_id in enumerate(items, 1):
        ent = syn.get(ent_id, downloadLocation=dest, ifcollision="keep.local")
        fh = ent._file_handle
        local = ent.path
        got = md5(local)
        if got != fh["contentMd5"]:
            raise SystemExit(f"md5 mismatch {ent_id} {local}: {got} vs {fh['contentMd5']}")
        rows.append((os.path.basename(local), f"synapse:{ent_id}.{ent.versionNumber}", str(os.path.getsize(local)),
                     got, "synapse_contentMd5:match", datetime.datetime.now().astimezone().isoformat(timespec="seconds")))
        print(f"{k}/{len(items)} {ent_id} {os.path.basename(local)}", flush=True)
    with open(os.path.join(dest, "manifest.tsv"), "w") as f:
        f.write("file\turl\tbytes\tmd5\tmd5_source\tretrieved\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"verified {len(rows)} / {len(items)}")


if __name__ == "__main__":
    main(sys.argv[1])
