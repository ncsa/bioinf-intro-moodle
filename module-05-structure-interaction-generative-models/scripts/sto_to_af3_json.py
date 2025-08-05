import sys
import json
import random

def parse_sto_to_af3_json(sto_path, output_json="af3_input.json"):
    seq_lines = {}
    ss_cons = ""

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

    pick_index = random.randint(0, N - 1)
    chosen_id = sequence_ids[pick_index]
    aligned_seq = seq_lines[chosen_id]



    print(f"Picking out the {pick_index}th sequence out of {N} sequences: {chosen_id}")

    # Remove gaps from both sequence and structure
    ungapped_seq = ""
    ungapped_ss = ""
    for s, b in zip(aligned_seq, ss_cons):
        if s != "-":
            ungapped_seq += s
            ungapped_ss += b

    # Convert dot-bracket to 0-based base pairs
    stack = []
    pairs = []
    for i, char in enumerate(ungapped_ss):
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            pairs.append([j, i])

    # Build AlphaFold3 input JSON
    af3_json = {
        "target": {
            "description": f"{chosen_id} from {sto_path}",
            "sequence": ungapped_seq,
            "molecule_type": "RNA"
        },
        "mode": "monomer",
        "template_mode": "none",
        "use_templates": False,
        "secondary_structure": pairs,
        "version": "1.0"
    }

    # Save with indentation for readability
    with open(output_json, 'w') as f:
	# Write all except secondary_structure nicely
        secondary_structure = af3_json.pop("secondary_structure")
        json.dump(af3_json, f, indent=2)
        f.write(",\n  \"secondary_structure\": [\n")
        compact = ", ".join(f"[{a},{b}]" for a, b in secondary_structure)
        f.write("    " + compact + "\n")
        f.write("  ]\n}\n")

    print(f"Saved AlphaFold3 input to {output_json}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sto_to_af3_json.py RFxxxxx.seed.sto")
    else:
        parse_sto_to_af3_json(sys.argv[1])

