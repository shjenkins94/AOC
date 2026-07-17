"""Rebuild a GARD nexus file with a larger alignment and new trees, using the original nexus file as a template."""

from Bio import AlignIO
from jinja2 import Template
import logging


def prepare_alignment_data(full_alignment_file):
    """Prepare alignment data from the full alignment file."""
    aln_data = {}

    logging.info(f"Preparing alignment data from: {full_alignment_file}")
    full_alignment = AlignIO.read(full_alignment_file, "fasta")

    tax_labels = [f"'{record.id}'" for record in full_alignment]
    tax_label_max_length = max(len(label) for label in tax_labels)

    aln_data["n_tax"] = len(tax_labels)
    logging.info(f"Number of taxa: {aln_data['n_tax']}")

    aln_data["tax_labels"] = " ".join(tax_labels)
    logging.info(f"Taxa labels: {aln_data['tax_labels']}")

    aln_data["n_char"] = full_alignment.get_alignment_length()
    logging.info(f"Alignment length (n_char): {aln_data['n_char']}")

    msa_matrix_lines = []
    for tax_label, record in zip(tax_labels, full_alignment):
        msa_matrix_lines.append(
            f"\t{tax_label.ljust(tax_label_max_length)} {record.seq}"
        )
    aln_data["msa_matrix"] = "\n".join(msa_matrix_lines)
    logging.info(f"New MSA matrix:\n{aln_data['msa_matrix']}")

    return aln_data


def extract_charset_lines(original_nexus_file):
    """Extract CHARSET lines from the original nexus file."""

    logging.info(f"Extracting CHARSET lines from: {original_nexus_file}")
    charset_lines = [line for line in open(original_nexus_file) if "CHARSET" in line]
    charsets = "".join(charset_lines)
    logging.info(f"CHARSET lines:\n{charsets}")

    return charsets


def format_tree_lines(tree_files):
    """Format tree lines from the tree files."""
    logging.info("Reading tree files")
    trees = []
    for tree_file in tree_files:
        with open(tree_file) as f:
            trees.append(f.read().strip())
    logging.info(f"Trees read: {len(trees)}")

    tree_lines = []
    for index, tree in enumerate(trees):
        tree_lines.append(f"\tTREE tree_{index + 1} = {tree}")

    tree_lines_str = "\n".join(tree_lines)
    logging.info(f"Tree lines:\n{tree_lines_str}")

    return tree_lines_str


def render_nexus_template(template_file, output_file, data):
    """Render the Jinja2 template and write to the output file."""
    logging.info(f"Rendering Jinja2 template: {template_file}")
    with open(template_file) as f:
        template = Template(f.read())

    logging.info(f"Writing output to: {output_file}")
    with open(output_file, "w") as f:
        f.write(template.render(**data))


def main():
    """Main Snakemake process."""
    logging.info("Starting main process")

    data = prepare_alignment_data(full_alignment_file=snakemake.input.full_alignment[0])

    data["charsets"] = extract_charset_lines(
        original_nexus_file=snakemake.input.bestgard[0]
    )

    data["trees"] = format_tree_lines(tree_files=snakemake.input.trees)

    render_nexus_template(
        template_file=snakemake.input.template,
        output_file=snakemake.output[0],
        data=data,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])
    main()
