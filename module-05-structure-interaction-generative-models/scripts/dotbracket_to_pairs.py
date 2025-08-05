import sys

def dotbracket_to_pairs(structure):
    stack = []
    pairs = []
    for i, char in enumerate(structure.strip()):
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            pairs.append([j, i])
    return pairs

def main(file_path):
    with open(file_path, 'r') as f:
        structure = f.read()
    pairs = dotbracket_to_pairs(structure)
    for pair in pairs:
        print(pair)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python dotbracket_to_pairs.py structure.txt")
    else:
        main(sys.argv[1])

