import searchsifter as ss
import os


def hash_location(o):
    location = os.path.join(o, "hashes.json")
    return location


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, default=100)
    parser.add_argument("-o", "--output-dir", type=str, default='')
    parser.add_argument("-p", "--pfam-filename", type=str)
    parser.add_argument("-t", "--pfam-file-type", type=str, choices=["regions", "stockholm"])
    args = parser.parse_args()
    if args.pfam_filename is not None:
        pfam = ss.relationships.pfam_file
        pfam_args = [args.pfam_filename]
        if args.pfam_file_type is not None:
            pfam_args.append(args.pfam_file_type)
    else:
        pfam = ss.relationships.pfam_db
        pfam_args = None
    hashes = ss.relationships.pfam.Hashes(n=args.n, pfam=pfam, pfam_args=pfam_args)
    hashes.save_to_file(hash_location(args.output_dir))
