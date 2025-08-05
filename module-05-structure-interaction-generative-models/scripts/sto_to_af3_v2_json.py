import sys
import json
import random

def parse_sto_to_af3_json(sto_path, output_json="../data/af3_input.json"):
    seq_lines = {}
    ss_cons = ""

    # Parse the .sto file
    with open(sto_path) as f:
        for line in f:
            if line.startswith("#=GC SS_cons"):
                ss_cons += line.strip().split()[2]
            elif line.startswith("#") or line.startswith("//"):
                continue
            elif line.strip():
                parts = line.strip().split()
                if len(parts) != 2:
                    continue
                name, seq = parts
                if name not in seq_lines:
                    seq_lines[name] = ""
                seq_lines[name] += seq

    sequence_ids = list(seq_lines.keys())
    N = len(sequence_ids)

    if N == 0:
        print("No sequences found.")
        return

    # Randomly pick a sequence
    pick_index = random.randint(0, N - 1)
    chosen_id = sequence_ids[pick_index]
    aligned_seq = seq_lines[chosen_id]

    print(f"Picking out the {pick_index}th sequence out of {N} sequences: {chosen_id}")

    # Remove alignment gaps
    ungapped_seq = ""
    ungapped_ss = ""
    for s, b in zip(aligned_seq, ss_cons):
        if s != "-":
            ungapped_seq += s
            ungapped_ss += b

    # Convert dot-bracket structure to pair restraints
    stack = []
    pair_restraints = []
    invalid = False
    for i, char in enumerate(ungapped_ss):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pair_restraints.append([j, i])
            else:
                invalid = True
                print(f"Warning: unmatched ')' at position {i}")
    if stack:
        print(f"Warning: {len(stack)} unmatched '(' remaining.")
        invalid = True

    # Build subunit
    subunit = {
        "sequence": ungapped_seq,
        "molecule_type": "RNA",
        "description": f"{chosen_id} from {sto_path}"
    }
    if not invalid and pair_restraints:
        subunit["pair_restraints"] = pair_restraints

    # Final AF3 JSON
    af3_json = {
        "fold_input": {
            "component_0": {
                "subunits": [subunit]
            }
        }
    }

    secondary = None
    if "pair_restraints" in subunit:
        secondary = subunit.pop("pair_restraints")

    with open(output_json, 'w') as f:
        json.dump(af3_json, f, indent=2)
        if secondary:
            f.seek(f.tell() - 31)  # go back just before "\n}"
            f.write(",\n          \"pair_restraints\": [")
            f.write(", ".join(f"[{a},{b}]" for a, b in secondary))
            f.write("]\n        }\n      ]\n    }\n  }\n}\n")

    print(f"Saved AlphaFold3 input to {output_json}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sto_to_af3_v2_json.py RFXXXXX.seed.sto")
    else:
        parse_sto_to_af3_json(sys.argv[1])

