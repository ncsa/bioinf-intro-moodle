import sys
import json

# Map each open bracket to its matching close bracket
bracket_pairs = {
    '(': ')',
    '[': ']',
    '{': '}',
    '<': '>',
    'A': 'a', 'B': 'b', 'C': 'c', 'D': 'd', 'E': 'e', 'F': 'f',
    'G': 'g', 'H': 'h', 'I': 'i', 'J': 'j', 'K': 'k', 'L': 'l',
    'M': 'm', 'N': 'n', 'O': 'o', 'P': 'p', 'Q': 'q', 'R': 'r',
    'S': 's', 'T': 't', 'U': 'u', 'V': 'v', 'W': 'w', 'X': 'x',
    'Y': 'y', 'Z': 'z'
}
# Create reverse map
closing_to_opening = {v: k for k, v in bracket_pairs.items()}

def dotbracket_to_pairs(structure):
    stacks = {k: [] for k in bracket_pairs}
    pairs = []
    for i, char in enumerate(structure.strip()):
        if char in bracket_pairs:
            stacks[char].append(i)
        elif char in closing_to_opening:
            open_char = closing_to_opening[char]
            if stacks[open_char]:
                j = stacks[open_char].pop()
                pairs.append([j, i])
    return sorted(pairs)

def main(file_path):
    with open(file_path, 'r') as f:
        structure = f.read()
    pairs = dotbracket_to_pairs(structure)
    print(json.dumps(pairs, indent=2))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python dotbracket_to_pairs.py structure.txt")
    else:
        main(sys.argv[1])

