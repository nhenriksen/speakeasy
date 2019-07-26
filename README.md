# speakeasy
Automates the application of SMIRNOFF99Frosst parameters to `mol2` files so they can be simulated with AMBER.

![Speakeasy logo.](speakeasy.png)




Make this:

```
speakeasy -i ligand.tripos.mol2 -fi -mol2 -o ligand.amber.mol2 -fo mol2 --smirnoff smirnoff99frosst.offxml --frcmod ligand.amber.frcmod
```

## Contributors

- Niel Henriksen (UCSD)
- David Slochower (UCSD)
- Jeff Wagner (Open Force Field Initiative)
- Jaime Rodríguez-Guerra (Charité - Universitätsmedizin Berlin)
- Rafal Wiewiora (MSKCC)
