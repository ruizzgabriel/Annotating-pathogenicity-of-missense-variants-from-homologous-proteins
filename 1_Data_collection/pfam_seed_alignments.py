# Extract individual Pfam seed alignments from the file 'Pfam-A.seed'. This file was retrieved directly from InterPro.


def extract_pfam_code(pfam_alignment):
    for line in pfam_alignment:
        if line.startswith('#=GF AC'):
            line = line.split(None, 3)
            pfam_code = line[-1].split('.')[-2]
            return pfam_code

def extract_sequences(pfam_alignment):
    sequences = []
    for line in pfam_alignment:
        if not line.startswith('#'):
            if '_HUMAN' in line:
                sequences.append(line)
    return sequences



pfam_seed = './data/Pfam-A.seed'
alignments = {}
seed_alignments = []

with open(pfam_seed, 'r', encoding='utf-8') as file:
    while True:
        try:
            # Save the content of the file
            line = next(file)
            if line.rstrip() != '//':
                seed_alignments.append(line)
            else:
                pfam_code = extract_pfam_code(seed_alignments)
                sequences = extract_sequences(seed_alignments)

                if sequences:
                    if pfam_code in alignments:
                        alignments[pfam_code] += sequences
                    else:
                        alignments[pfam_code] = sequences
                    #print(pfam_code)
                    #print(sequences)

                seed_alignments = []
        except UnicodeDecodeError:
            continue
        except StopIteration:
            break

#print(len(alignments))

output_path = './data/Pfam_SEED_HUMAN_2024_03_20'
for key in alignments:
    with open(f'{output_path}/{key}_HUMAN.txt', 'w+') as output_file:
        output_file.write(''.join(alignments[key]))