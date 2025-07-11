# cnn_sanity.py
def create_tsv(filename, data, inputs, outputs):
    with open(filename, 'w') as f:
        f.write(f"#inputs\t{inputs}\n")
        f.write(f"#outputs\t{outputs}\n")
        for row in data:
            f.write('\t'.join(map(str, row)) + '\n')

# Simple sanity check
sanity_data = [
    [1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
]
create_tsv('cnn_sanity.tsv', sanity_data, 5, 1)

# XOR test
xor_data = [
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 0.0]
]
create_tsv('cnn_xor.tsv', xor_data, 2, 1)