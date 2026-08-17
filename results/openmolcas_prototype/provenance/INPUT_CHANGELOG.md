# Dy input change record

Parent input: official OpenMolcas `test/additional/321.input` at commit `8355057f32d65706a35996b5ab07cac2962bb728`.

Validated prototype changes:

1. Replaced `NoCD` with explicit `RICD` and `CDTH 1.0D-6`.
2. Removed the legacy explicit `HINT 0.0 5.0 21`, `TMAG 1 2.0`, and `TYPE 7` lines. The core MLTP/MVEC SINGLE_ANISO request remained, and OpenMolcas used its printed default magnetic-property grid.
3. Added no project geometry, no manuscript value, and no estimated scientific data.

The original official input and the exact validated derivative are both in `inputs/`.
