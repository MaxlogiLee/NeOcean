from Bio import SeqIO

def reverse_fasta_sequences(input_file, output_file):
  with open(input_file, "r") as in_handle:
    with open(output_file, "w") as out_handle:
      for record in SeqIO.parse(in_handle, "fasta"):
        reversed_seq = record.seq[::-1]
        reversed_record = record
        reversed_record.seq = reversed_seq
        SeqIO.write(reversed_record, out_handle, "fasta")

input_file="uniprot_merge_protein.fa"
output_file="Decoy_uniprotkb.fa"
reverse_fasta_sequences(input_file, output_file)
