import os
from searchsifter import relationships
import json


def exact_jaccard(test_accs, pfam_accs, family_source):
    for pfam_acc in pfam_accs:
        pfam_fam = family_source(pfam_acc)
        for test_acc in test_accs:
            test_fam = family_source(test_acc)
            intersect = sum(len(r) for _, rs in pfam_fam.overlap(test_fam).items() for r in rs)
            union = sum(len(r) for _, rs in pfam_fam.union(test_fam).items() for r in rs)
            ji = intersect / union
            jc = intersect / pfam_fam.residues_covered()
            yield pfam_acc, test_acc, ji, jc


def est_jaccard(test_accs, hash_paths, ns, family_source):
    for path in hash_paths:
        dir_, name = os.path.split(path)
        window = int(name.split('_')[1].split('.')[0])
        for n in ns:
            with open(path) as hash_file:
                hashes = relationships.pfam.ResidueHashes(window, n=n, from_file=hash_file)
            for test_acc in test_accs:
                fam = family_source(test_acc)
                jis = hashes.estimate_jaccard(fam)
                jcs = hashes.estimate_containment(fam)
                for pfam_acc, ji in jis.items():
                    jc = jcs[pfam_acc]
                    yield pfam_acc, test_acc, ji, jc, n, window


if __name__ == "__main__":
    import argparse
    from glob import glob

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--ns", nargs='+', type=int)
    parser.add_argument("-a", "--hashes", nargs='+', type=str)
    parser.add_argument("-s", "--samples", type=str)
    parser.add_argument("-p", "--pfam-filename", type=str)
    parser.add_argument("-t", "--pfam-file-type", type=str, choices=["regions", "stockholm"])

    args = parser.parse_args()

    with open(args.samples) as sample_file:
        sample = json.load(sample_file)

    hash_paths = [p for p in args.hashes for q in glob(p)]

    with open(hash_paths[0]) as hash_file:
        all_pfam = json.load(hash_file).keys()

    if args.pfam_filename is not None:
        pfam_db = relationships.pfam_file.PfamFamilies(args.pfam_filename, args.pfam_file_type, all_pfam)
        family_source = pfam_db.__getitem__
    else:
        family_source = relationships.pfam_db.PfamFamily.from_accession

    if args.ns is None:
        for cols in exact_jaccard(sample, all_pfam, family_source):
            print(*cols, 1, 1, "exact", sep='\t')
    else:
        for cols in est_jaccard(sample, hash_paths, args.ns, family_source):
            print(*cols, "estimated", sep='\t')
