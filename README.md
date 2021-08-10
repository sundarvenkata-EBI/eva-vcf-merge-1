# EVA VCF Merge

Library to merge VCF files horizontally and vertically.

## Horizontal Merge

Horizontally mergeable VCFs have non-overlapping sample sets and can be combined to create a single multi-sample
file, using e.g. [`bcftools merge`](http://samtools.github.io/bcftools/bcftools.html#merge).

## Vertical Merge

Vertically mergeable VCFs have identical sample sets and can be combined to create a single multi-chromosome or
multi-variant file, using e.g. [`bcftools concat`](http://samtools.github.io/bcftools/bcftools.html#concat).
