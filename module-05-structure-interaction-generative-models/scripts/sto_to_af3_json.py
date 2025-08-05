import sys
import json
import random

def parse_dot_bracket(dot_bracket, verbose=True):
    stack = []
    pairs = []
    for i, c in enumerate(dot_bracket.strip()):
        if c in "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            stack.append((c, i))
        elif c in ")]}>abcdefghijklmnopqrstuvwxyz":
            if not stack:
                if verbose:
                    print(f"Warning: unmatched '{c}' at position {i}")
                continue
            open_char, j = stack.pop()
            pairs.append([j, i])
    return pairs

class CompactArraysEncoder(json.JSONEncoder):
    def default(self, o):
        return super().default(o)

    def encode(self, obj):
        json_str = super().encode(obj)

        # Flatten only the "pair_restraints" array of arrays
        import re
        pattern = r'"pair_restraints": \[\s*(\[[^\[\]]+\](?:,\s*\[[^\[\]]+\])*)\s*\]'
        def replacer(match):
            flat = match.group(1).replace('\n', '').replace('  ', '')
            return f'"pair_restraints": [\n  {flat}\n]'

        return re.sub(pattern, replacer, json_str)


def parse_sto_to_af3_json(sto_file, output_json="../data/rfam_example.json"):
    with open(sto_file) as f:
        lines = f.readlines()

    sequence_lines = []
    ss_cons_line = ""
    for line in lines:
        if not line.startswith("#") and line.strip():
            sequence_lines.append(line.strip())
        elif line.startswith("#=GC SS_cons"):
            ss_cons_line = line.strip().split(maxsplit=2)[-1]

    if not sequence_lines:
        raise ValueError("No sequences found.")

    N = len(sequence_lines)
    pick_index = random.randint(0, N - 1)
    picked_line = sequence_lines[pick_index]
    seq_id, seq = picked_line.split(maxsplit=1)

    print(f"Picking out the {pick_index}th sequence out of {N} sequences: {seq_id}")

    subunit = {
        "sequence": seq.replace('.', '').replace('-', ''),
        "molecule_type": "RNA",
        "description": f"{seq_id} from {sto_file}"
    }

    if ss_cons_line:
        pair_restraints = parse_dot_bracket(ss_cons_line)
        if pair_restraints:
            subunit["pair_restraints"] = pair_restraints

    af3_json = {
        "component_0": {
            "subunits": [subunit]
        }
    }

    with open(output_json, 'w') as f:
        f.write(json.dumps(af3_json, indent=2, cls=CompactArraysEncoder))

if __name__ == "__main__":
    parse_sto_to_af3_json(sys.argv[1])

