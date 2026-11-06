# Input files

The script reads these files from this folder, using these exact names. They were not archived when the paper was published, so [ADD: either the recovered original files, or a note that they were reconstructed from GenBank sequences and the paper's methods, and may differ slightly from the published analysis].

| File | Content |
|---|---|
| `tree.nwk` | Time-calibrated phylogeny in Newick format, branch lengths in millions of years |
| `geography.data` | Tip ranges in Lagrange PHYLIP format |
| `timeperiods.txt` | Time slices (ages in millions of years) |
| `manual_dispersal_multipliers.txt` | Dispersal multipliers per time slice |
| `areas_allowed.txt` | Areas allowed in each time slice |
| `areas_adjacency.txt` | Area adjacency per time slice |
| `distances_matrix.txt` | Distances between areas |

## Format notes

**Phylogeny:** binary (no polytomies), no negative branch lengths, no spaces in tip names (use `_`). Very short branches are treated as direct ancestors.

**Geography file:** first line = number of tips, a tab, the number of areas and the area names in parentheses. Then one line per tip: name, a tab, and the 0/1 presence string with no spaces.

```
10	4 (A B C D)
species_A	0110
species_B	0111
species_C	0001
```

Tip names must match the phylogeny exactly. Tips should be genetically distinct lineages, not individual specimens.
