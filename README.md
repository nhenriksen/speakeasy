# speakeasy
Automates the conversion of SMIRNOFF parameters to AMBER force field files

Make this:

```
speakeasy -i ligand.tripos.mol2 -fi -mol2 -o ligand.amber.mol2 -fo mol2 --smirnoff smirnoff99frosst.offxml --frcmod ligand.amber.frcmod
```