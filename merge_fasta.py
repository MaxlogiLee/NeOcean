def preprocess_fasta(fasta_path, temp_fasta_path):
    with open(fasta_path, 'r') as infile, open(temp_fasta_path, 'w') as outfile:
        sequence = ""
        for line in infile:
            if line.startswith(">"):
                if sequence:
                    outfile.write(sequence + "\n")
                    sequence = ""
                outfile.write(line)  
            else:
                sequence += line.strip()
        if sequence:
            outfile.write(sequence + "\n")

if __name__ == "__main__":
    input_file = "sample.fasta"
    
    output_file = "sample_merge.fasta"
    preprocess_fasta(input_file, output_file)
    print(f"Processing complete! Results have been saved to: {output_file}")
