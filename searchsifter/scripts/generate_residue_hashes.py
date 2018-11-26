import searchsifter.relationships.pfam as pf
from searchsifter import relationships
import os
import time


def hash_location(o, w):
    return os.path.join(o, "rhashes_{}.json".format(w))


def size_location(o, w):
    return os.path.join(o, "sizes_{}.json".format(w))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, required=True)
    parser.add_argument("-w", "--windows", nargs='+', type=int, required=True)
    parser.add_argument("-o", "--output-dir", type=str, default='')
    parser.add_argument("-p", "--pfam-filename", type=str)
    parser.add_argument("-t", "--pfam-file-type", type=str)
    args = parser.parse_args()
    if args.pfam_filename is not None:
        pfam = relationships.pfam_file
        pfam_args = [args.pfam_filename]
        if args.pfam_file_type is not None:
            pfam_args.append(args.pfam_file_type)
    else:
        pfam = relationships.pfam_db
        pfam_args = None
    hashes = pf.ResidueHashes.hashes_with_windows(args.windows, args.n, pfam=pfam, pfam_args=pfam_args)
    for w, h in hashes.items():
        t = time.time()
        h.save_to_file(hash_location(args.output_dir, w))
        print("Generated hash with w={} in {} seconds".format(w, time.time() - t))
