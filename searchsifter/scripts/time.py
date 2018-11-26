from searchsifter import relationships
import json
import os
import timeit


def inner_time_minhash(test_accs, hash_path, n, window, family_source):
    with open(hash_path) as hash_file:
        hashes = relationships.pfam.ResidueHashes(window, n=n, from_file=hash_file)
        for test_acc in test_accs:
            fam = family_source(test_acc)
            ji_time = min(timeit.repeat("hashes.estimate_jaccard(fam)",
                                        globals=locals(),
                                        number=1,
                                        repeat=3))
            jc_time = min(timeit.repeat("hashes.estimate_containment(fam)",
                                        globals=locals(),
                                        number=1,
                                        repeat=3))
            yield test_acc, ji_time, jc_time, n, window, len(fam.proteins())


def time_minhash(test_accs, hash_paths, ns, family_source):
    for path in hash_paths:
        dir_, name = os.path.split(path)
        window = int(name.split('_')[1].split('.')[0])
        for n in ns:
            yield from inner_time_minhash(test_accs, path, n, window, family_source)


if __name__ == '__main__':
    import argparse
    from glob import glob

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--ns", nargs='+', type=int)
    parser.add_argument("-a", "--hashes", nargs='+', type=str)
    parser.add_argument("-s", "--samples", type=str)
    parser.add_argument("-p", "--pfam-filename", type=str)
    parser.add_argument("-t", "--pfam-file-type", type=str, choices=["regions", "stockholm"])

    args = parser.parse_args()

    if args.pfam_filename is not None:
        pfam_db = relationships.pfam_file.PfamFamilies(args.pfam_filename, args.pfam_file_type)
        family_source = pfam_db.__getitem__
    else:
        family_source = relationships.pfam_db.PfamFamily.from_accession

    with open(args.samples) as sample_file:
        sample = json.load(sample_file)

    hash_paths = [p for p in args.hashes for q in glob(p)]

    for cols in time_minhash(sample, hash_paths, args.ns, family_source):
        print(*cols, sep='\t')
