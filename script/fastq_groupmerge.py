#!/usr/bin/env python
import gzip
import shutil
import pandas as pd
from pathlib import Path
import argparse
from script.version import __version__


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Merge paired or single FASTQ (including FASTQ-Sanger) files based on metadata, samples can belong to multiple groups, or merge all reads together into one file."
    )
    parser.add_argument(
        "--metadata",
        help="Path to metadata CSV/TSV file. If no metadata is included all files in fastq_dir will be merged to one forward and one reverse read!",
    )
    parser.add_argument("fastq_dir", help="Directory containing FASTQ files")
    parser.add_argument("output_dir", help="Output directory for merged FASTQs")
    parser.add_argument(
        "--group_col", default="group", help="Metadata column name for grouping"
    )
    parser.add_argument(
        "--sep", default=",", help="Column separator in metadata (default: ',')"
    )
    parser.add_argument(
        "--forward_suffix",
        default="_forward",
        help="Suffix to find the forward reads (default: _forward)",
    )
    parser.add_argument(
        "--reverse_suffix",
        default="_reverse",
        help="Suffix to find the reverse reads (default: _reverse)",
    )
    parser.add_argument("-v", "-V", "--version", action="version", version=__version__)
    parser.add_argument(
        "--single_reads",
        action="store_true",
        help="Set this flag if only single reads are inputted for this tool to merge!",
    )
    parser.print_usage = parser.print_help

    args = parser.parse_args()
    return args


def find_fastq_files(in_dir, file_name):
    """Return possible FASTQ/FastqSanger files with different extensions."""
    # ✅ Include all FASTQ variants
    extensions = [".fastq.gz", ".fastqsanger.gz", ".fastq", ".fastqsanger"]
    for ext in extensions:
        file = in_dir / f"{file_name}{ext}"
        if file.exists():
            return file
    return None


def merge_fastqs_pair(
    metadata_file,
    fastq_dir,
    output_dir,
    group_col="group",
    sep=",",
    forward_suffix="_forward",
    reverse_suffix="_reverse",
):
    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load metadata
    df = pd.read_csv(metadata_file, sep=sep)
    if "sample_id" not in df.columns or group_col not in df.columns:
        raise ValueError(f"Metadata must have columns: sample_id and {group_col}")

    df.dropna(inplace=True)

    groups = df.groupby(group_col)["sample_id"].apply(list).to_dict()

    for group, samples in groups.items():
        print(f"\n🔹 Merging samples for group '{group}' -> {samples}")

        out_R1 = output_dir / f"{group}{forward_suffix}.fastq.gz"
        out_R2 = output_dir / f"{group}{reverse_suffix}.fastq.gz"

        for out_file in (out_R1, out_R2):
            if out_file.exists():
                out_file.unlink()

        for sample in samples:
            R1 = find_fastq_files(fastq_dir, f"{sample}{forward_suffix}")
            R2 = find_fastq_files(fastq_dir, f"{sample}{reverse_suffix}")

            if not R1 or not R2:
                print(f"⚠️  Skipping {sample}: missing one of the paired files.")
                continue

            print(f"  ➕ Adding {R1.name} and {R2.name} to group {group}")
            for src, dest in [(R1, out_R1), (R2, out_R2)]:
                open_func = gzip.open if src.suffix.endswith("gz") else open
                with open_func(src, "rb") as f_in, gzip.open(dest, "ab") as f_out:
                    shutil.copyfileobj(f_in, f_out)

        print(f"✅ Done: {out_R1.name}, {out_R2.name}")

    print("\n🎉 All merges complete!")


def merge_all_pair(fastq_dir, output_dir, forward_suffix="_forward", reverse_suffix="_reverse"):
    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # ✅ Include .fastqsanger variants too
    file_list = (
        list(fastq_dir.glob("*.fastq.gz"))
        + list(fastq_dir.glob("*.fastqsanger.gz"))
        + list(fastq_dir.glob("*.fastq"))
        + list(fastq_dir.glob("*.fastqsanger"))
    )

    if not file_list:
        print("❌ No FASTQ or FASTQ-Sanger files found.")
        return

    forward_reads = sorted([f for f in file_list if forward_suffix in f.name])
    reverse_reads = sorted([f for f in file_list if reverse_suffix in f.name])

    print(f"📂 Found {len(forward_reads)} forward and {len(reverse_reads)} reverse files.")

    out_R1 = output_dir / f"merged{forward_suffix}.fastq.gz"
    out_R2 = output_dir / f"merged{reverse_suffix}.fastq.gz"

    for out_file in (out_R1, out_R2):
        if out_file.exists():
            out_file.unlink()

    print(f"🧬 Merging {len(forward_reads)} files into {out_R1.name}")
    for read in forward_reads:
        print(f"➕ Adding {read.name}")
        open_func = gzip.open if read.suffix.endswith("gz") else open
        with open_func(read, "rb") as f_in, gzip.open(out_R1, "ab") as f_out:
            shutil.copyfileobj(f_in, f_out)

    print(f"🧬 Merging {len(reverse_reads)} files into {out_R2.name}")
    for read in reverse_reads:
        print(f"➕ Adding {read.name}")
        open_func = gzip.open if read.suffix.endswith("gz") else open
        with open_func(read, "rb") as f_in, gzip.open(out_R2, "ab") as f_out:
            shutil.copyfileobj(f_in, f_out)

    print("\n🎉 All merges complete!")


def merge_fastqs_single(metadata_file, fastq_dir, output_dir, group_col="group", sep=","):
    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    df = pd.read_csv(metadata_file, sep=sep)
    if "sample_id" not in df.columns or group_col not in df.columns:
        raise ValueError(f"Metadata must have columns: sample_id and {group_col}")

    df.dropna(inplace=True)

    groups = df.groupby(group_col)["sample_id"].apply(list).to_dict()

    for group, samples in groups.items():
        print(f"\n🔹 Merging samples for group '{group}' -> {samples}")

        out = output_dir / f"{group}.fastq.gz"
        if out.exists():
            out.unlink()

        for sample in samples:
            R = find_fastq_files(fastq_dir, sample)
            if not R:
                print(f"⚠️  Skipping {sample}: missing file.")
                continue

            print(f"  ➕ Adding {R.name} to group {group}")
            with (gzip.open if R.suffix.endswith("gz") else open)(R, "rb") as f_in, gzip.open(out, "ab") as f_out:
                shutil.copyfileobj(f_in, f_out)

        print(f"✅ Done: {out.name}")

    print("\n🎉 All merges complete!")


def merge_all_single(fastq_dir, output_dir):
    fastq_dir = Path(fastq_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # ✅ Include .fastqsanger variants too
    file_list = (
        list(fastq_dir.glob("*.fastq.gz"))
        + list(fastq_dir.glob("*.fastqsanger.gz"))
        + list(fastq_dir.glob("*.fastq"))
        + list(fastq_dir.glob("*.fastqsanger"))
    )

    if not file_list:
        print("❌ No FASTQ or FASTQ-Sanger files found.")
        return

    out = output_dir / "merged.fastq.gz"
    if out.exists():
        out.unlink()

    for read in file_list:
        print(f"➕ Adding {read.name}")
        open_func = gzip.open if read.suffix.endswith("gz") else open
        with open_func(read, "rb") as f_in, gzip.open(out, "ab") as f_out:
            shutil.copyfileobj(f_in, f_out)

    print("\n🎉 All merges complete!")


if __name__ == "__main__":
    args = parse_arguments()

    if args.single_reads:
        if args.metadata:
            merge_fastqs_single(
                args.metadata,
                args.fastq_dir,
                args.output_dir,
                group_col=args.group_col,
                sep=args.sep,
            )
        else:
            merge_all_single(args.fastq_dir, args.output_dir)
    else:
        if args.metadata:
            merge_fastqs_pair(
                args.metadata,
                args.fastq_dir,
                args.output_dir,
                group_col=args.group_col,
                sep=args.sep,
                forward_suffix=args.forward_suffix,
                reverse_suffix=args.reverse_suffix,
            )
        else:
            merge_all_pair(
                args.fastq_dir,
                args.output_dir,
                forward_suffix=args.forward_suffix,
                reverse_suffix=args.reverse_suffix,
            )
