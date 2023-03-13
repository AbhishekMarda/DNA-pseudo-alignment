from Bio import SeqIO
import pandas as pd


K_VALUE = 21
ISOFORM_READ_FILE = "chr11_transcriptome.fasta"
TRANSCRIPT_READ_FILE = "reads.fasta"


isoform_file = open(ISOFORM_READ_FILE, 'r')

# dictionary indexed by a k_mer pointing back to a set of isoforms it could have come from
back_index = dict()

for isoform in SeqIO.parse(isoform_file, "fasta"):
    seq = isoform.seq
    seq_str = str(seq)
    seq_len = len(seq_str) 
    for i in range(0, seq_len - K_VALUE + 1):
        k_mer = seq_str[i:i+K_VALUE]
        if k_mer in back_index:
            back_index[k_mer].add(isoform.id)
        else:
            back_index[k_mer] = {isoform.id}
    
    
isoform_file.close()

# create a function such that given the reverse maps, we can find 
# the equivalence class for the transcript
def get_equivalence_class(isoform_dict, transcript) -> set | None:
    transcript_str = str(transcript) 
    transcript_len = len(transcript_str)
    
    equivalence_class = None
    # take intersection over all the equivalence sets for each k-mer
    for i in range(0, transcript_len - K_VALUE + 1):
        k_mer = transcript_str[i:i+K_VALUE]

        # check if the k_mer exists in the dictionary
        if k_mer in isoform_dict:
            isoforms = isoform_dict[k_mer]
        else:
            isoforms = set()

        if not equivalence_class: # first iteration, no set to intersect on yet
            equivalence_class = isoforms
        else: 
            equivalence_class = equivalence_class.intersection(isoforms)
    return equivalence_class


transcript_file = open(TRANSCRIPT_READ_FILE, 'r')

equivalence_class_mappings = dict()

for transcript in SeqIO.parse(transcript_file, "fasta"):
    seq = transcript.seq
    # first work with the forward strand
    canonical_strand_classes = get_equivalence_class(back_index, seq)
    if canonical_strand_classes is None: 
        break 
    # get the classes for the reverse strand
    reverse_complement_classes = get_equivalence_class(back_index, seq.reverse_complement())
    if reverse_complement_classes is None:
        break
    """
    Current implementation for difference in forward and reverse
    strands is to take the union
    """
    combined_classes = canonical_strand_classes.union(reverse_complement_classes)
    
    hashable_combined_classes = frozenset(combined_classes)
    if hashable_combined_classes in equivalence_class_mappings:
        equivalence_class_mappings[hashable_combined_classes] += 1
    else:
        equivalence_class_mappings[hashable_combined_classes] = 1
transcript_file.close()

# sort by length of the equivalent sets
temp_dict_list = sorted(list(equivalence_class_mappings.items()), key = lambda key : len(key[0]))

equivalence_class_mappings = {ele[0] : ele[1] for ele in temp_dict_list}

data = list()
for equivalence_class, count in equivalence_class_mappings.items():
    data.append([count, len(equivalence_class), list(equivalence_class)])


output_df = pd.DataFrame(data, columns=["counts", "num items in equivalence classe", "isoforms in equivalence class"])

with open("output.csv", 'w') as output_file:
    output_file.write(output_df.to_csv(index=False))