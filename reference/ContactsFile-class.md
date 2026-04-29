# `ContactsFile` S4 class

The `ContactsFile` class describes a `BiocFile` object, pointing to the
location of an Hi-C matrix file (cool, mcool, hic, hicpro, ...) and
containing additional slots:

1.  resolution: at which resolution the associated mcool file should be
    parsed

2.  pairsFile: the path (in plain character) to an optional pairs file
    (stored as a `PairsFile` object);

3.  metadata: a list. If the CoolFile is created by `HiCool`, it will
    contain two elements: `log` (path to `HiCool` processing log file)
    and `stats` (aggregating some stats from `HiCool` mapping).

ContactsFile methods.

## Arguments

- path:

  String; path to an Hi-C matrix file (cool, mcool, hic, hicpro)

- resolution:

  numeric; resolution to use with Hi-C matrix file

- pairsFile:

  String; path to a pairs file

- metadata:

  list.

- object:

  A `ContactsFile` object.

- x:

  A `ContactsFile` object.

## Slots

- `resolution`:

  numeric value or NULL

- `pairsFile`:

  PairsFile object

- `metadata`:

  list

## See also

[`CoolFile()`](CoolFile-class.md), [`HicFile()`](HicFile-class.md),
[`HicproFile()`](HicproFile-class.md)
