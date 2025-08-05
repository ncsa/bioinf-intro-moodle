## character lengths to subtract when reformating the pair_restraints line

in code developing, we noticed that when we directly dump pair_restraints into the json file, it will look like:
```bash
"pair_restraints": [
  [
    0,
    71
  ],
  [
    1,
    70
  ]
]
```
Therefore, we pop it from the end of the subunit, and add it back after we've dumped other parts into the json file.

In order to do this, we need to remove the ending "]"s and "}"s, add the formatted pair restraints, and pad the "]"s and "}"s back. We need to calculate how much characters to rewind.

```bash
len("        }\n") = 9
len("      ]\n") = 8
len("    }\n") = 7
len("  }\n") = 5
len("}\n") = 2
-----------------
Total = 31 characters
```

Therefore, we use f.seek(f.tell()-31) to find the place to insert the pair_restraints. 
