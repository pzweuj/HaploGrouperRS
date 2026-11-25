deCODE fork of ISOGG 2019 14.196 (isogg2019-decode1)

With haploGrouper v1.1, we released a forked version of the ISOGG 2019 14.196
Y-chromosome locus file. We recommend using this isogg2019-decode1 locus file
instead of the original isogg2019 locus file as it can lead to substantially
more detailed haplogroup assignments in some use cases.

It is important to note that we made no changes to the tree file, so the
original ISOGG 14.196 phylogeny structure and haplogroup names are still used
for assignment.

We made our changes after inspecting the mismatchLoci detected in the
assignments of around 300,000 Y chromosomes available at deCODE genetics. This
revealed many variants - more precisely, mutant alleles resulting from
mutations at specified branches - with a high mismatch rate across samples,
suggesting that they did not truly belong to the branch they were assigned to
in ISOGG. 

Depending on how HaploGrouper was run, this could result in less finely resolved
haplogroup assignments for samples than expected: some branches carried a large
number of these frequently mismatching variants, and individuals would often
get "stuck" behind the mismatchScore burden of that branch and be prevented
from progressing unless they could rack up a compensating amount of supporting
variants downstream.

For each variant with a high mismatch rate, we evaluated whether there was a
downstream branch where we could place the mutation event to result in a <1%
mismatch rate. If there wasn’t, we discarded the variant. We then reassigned
haplogroups using the edited locus file. This procedure was repeated once to
allow previously blocked-off haplogroups (as described above) to be assessed.
Of the original 71405 mutations listed in the file, we discarded 3069 and moved
1170 farther downstream. We were unable to evaluate 3003 mutations because no
sample was assigned to associated branches. 

Branches with many high-mismatch rate variants were often for haplogroups ending
with a tilde symbol (~). Around two-thirds of all the variants we moved or
discarded were on ~ branches. In ISOGG, this symbol represents some kind of
ambiguity, and seems to have been employed in different ways by different ISOGG
volunteers. Perhaps the branching around that node was uncertain, or such
haplogroups were used to collect mutations believed to belong on some as-yet
unascertained position downstream. When we first compiled the locus file for
ISOGG 2019, we used the tree as-is and gave no special treatment to ~
variants. 

ISOGG also has uncertainty indicators for individual variants, "^" and "^^".
They are typically not applied to variants on ~ branches. After excluding ~
variants, ^/^^ variants are more likely to get removed: 8.5% of ^^, 4% of ^,
and 3% of unmarked variants were discarded. There are also minor differences in
chance to be moved downstream: 1.5% of ^^, 4% of ^, and 2% of unmarked
variants.

The impact of these high-mismatch rate variants depends on the properties of the
genotype data used as input. Firstly, many of these variants were very rare and
would typically not appear in VCFs that only include polymorphic loci.
Secondly, some haplogroups were more affected than others, meaning that the
ancestry of the samples analysed mediates impact. For example, many samples
were blocked from progressing to and beyond R1b1a1b1a1a2, a haplogroup very
common in western Europe and especially Ireland (where it shows ~60%
frequency).

We can illustrate this disparate impact with two examples. When using the new
locus file to re-assign haplogroups for the 1000G VCF, which only includes loci
that are polymorphic in that set, and without using the -f reference sequence
parameter, 6.7% of assignments change, and the average assignment is 0.1
branches deeper (more finely resolved). However, in our analysis of UKB samples
using all sequenced loci, 45% of assignments change, and the average assignment
is 2.9 branches deeper.

Our version of the tree is forked from ISOGG 14.196, which is not the most
recent (and likely final ever) version, 15.73. This means that our tree is
missing some subsequent ISOGG changes. These mainly affect G and C, but some
edits were also made to D, M, S, and P. Our Y-chromosome data is not an even
sampling of the tree – many R's, some A's, just 20 or so M's – so some
haplogroups were better evaluated than others. However, our biases are probably
similar to ISOGG's, so the sections we check less well were probably less
elaborated in the original trees.

Due to our changes, some branches no longer have any mutation events placed on
them. This isn't a huge problem because haploGrouper evaluates support for
every haplogroup in the tree, so samples can still be assigned downstream from
a branch with no mutation events.