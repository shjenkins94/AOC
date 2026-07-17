from Bio import AlignIO
from jinja2 import Template
import logging

# Input files
# full_alignment_file = "tiny.SA.codons.cln.fa"
# tree_files = ["part1/FastTree.labelled.treefile", "part2/FastTree.labelled.treefile", "part3/FastTree.labelled.treefile", "part4/FastTree.labelled.treefile"]
# original_nexus_file = "tiny.RD.SA.codons.cln.fa.cluster.fasta.best-gard"
# template_file = "best-gard.jinja2"
# output_file = "best-gard.nexus"

full_alignment_file = snakemake.input["full_alignment"]
tree_files = snakemake.input["trees"]
original_nexus_file = snakemake.input["bestgard"]
template_file = snakemake.input["template"]
output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Set up logging
logging.basicConfig(level=logging.INFO, filename=log_file)

# read full alignment file to get n_tax, tax_labels, n_char, and msa_matrix
logging.info(f"Reading full alignment file: {full_alignment_file}")
full_alignment = AlignIO.read(full_alignment_file, "fasta")

tax_labels = [f"'{record.id}'" for record in full_alignment]
tax_label_max_length = max(len(label) for label in tax_labels)
n_char = full_alignment.get_alignment_length()
logging.info(f"Number of taxa: {len(tax_labels)}")
logging.info(f"Taxa labels: {tax_labels}")
logging.info(f"Alignment length (n_char): {n_char}")

logging.info("Building MSA matrix lines")
msa_matrix_lines = []
for tax_label, record in zip(tax_labels, full_alignment):
    msa_matrix_lines.append(f"\t{tax_label.ljust(tax_label_max_length)} {record.seq}")

logging.info(f"New MSA matrix:")
logging.info('\n'.join(msa_matrix_lines))

logging.info(f"Reading CHARSET lines from original nexus file: {original_nexus_file}")
# Read the CHARSET lines from the original nexus file
charsets = [line for line in open(original_nexus_file) if "CHARSET" in line]

logging.info(f"CHARSET lines:")
logging.info(''.join(charsets))

# Read each tree file to get trees
logging.info("Reading tree files")
trees = []
for tree_file in tree_files:
    with open(tree_file) as f:
        trees.append(f.read().strip())
logging.info(f"Trees read: {len(trees)}")

tree_lines = []
for index, tree in enumerate(trees):
    tree_lines.append(f"\tTREE tree_{index + 1} = {tree}")

logging.info(f"Tree lines:")
logging.info('\n'.join(tree_lines))

# Prepare data for the Jinja2 template
data = {
    "n_tax": len(tax_labels),
    "tax_labels": ' '.join(tax_labels),
    "n_char": n_char,
    "msa_matrix": '\n'.join(msa_matrix_lines),
    "charsets": ''.join(charsets),
    "trees": '\n'.join(tree_lines)
}

# Render the Jinja2 template and write to the output file
logging.info(f"Rendering Jinja2 template: {template_file}")
with open(template_file) as f:
    template = Template(f.read())

logging.info(f"Writing output to: {output_file}")
with open(output_file, "w") as f:
    f.write(template.render(**data))