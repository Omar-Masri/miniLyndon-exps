import sys
import statistics
import csv
import io


def fasta_lengths(stream):
    seq_len = 0
    for line in stream:
        if line.startswith(">"):
            if seq_len > 0:
                yield seq_len
                seq_len = 0
        else:
            seq_len += len(line.strip())
    if seq_len > 0:
        yield seq_len


def fastq_lengths(stream):
    i = 0
    for line in stream:
        i += 1
        if i % 4 == 2:
            yield len(line.strip())


def main():
    first = sys.stdin.readline()
    if first.startswith(">"):
        lengths = list(fasta_lengths([first] + sys.stdin.readlines()))
    elif first.startswith("@"):
        lengths = list(fastq_lengths([first] + sys.stdin.readlines()))
    else:
        print("Input is not in FASTA or FASTQ format", file=sys.stderr)
        return

    if not lengths:
        print("No sequences found.")
        return
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(
        [
            "min",
            "max",
            "median",
            "mean",
            "stdev",
            "total_sequences",
            "total_bases",
        ]
    )
    writer.writerow(
        [
            min(lengths),
            max(lengths),
            statistics.median(lengths),
            f"{sum(lengths)/len(lengths):.2f}",
            f"{statistics.stdev(lengths):.2f}",
            len(lengths),
            sum(lengths),
        ]
    )
    print(output.getvalue(), end="")


if __name__ == "__main__":
    main()
