# IGV practice: validating your variants

By now you have called variants, filtered them, annotated them, and picked out a candidate mutation for FMF. Everything so far has been a table — a caller's opinion about what is in your data, summarised into a line of text.

In this session we go back and look at the **reads themselves**. This is the last check before you believe a result: a variant call is only as good as the alignment underneath it, and some things — a variant supported by two reads, a region with no coverage at all, a call that only appears because of a misalignment — are invisible in a VCF but obvious the moment you look.

---

## Input Data

IGV runs on your own computer, so the files have to come down from the server. You need **six** of them — two BAMs with their indexes, and one VCF with its index:

| File | Where it comes from |
|---|---|
| `wes46_dedup.bam` (+ its index) | `data/bam_clean/` — the alignment of sample 1 |
| `wes78_dedup.bam` (+ its index) | `data/bam_clean/` — the alignment of sample 2 |
| `filtered_snps.rsid.vcf.gz` (+ `.tbi`) | `data/annotation/` — your annotated variants, with rsIDs |

Both BAMs came out of [the calling pipeline](../varinat_calling/variant_calling_pipeline.md) (Step 3, `MarkDuplicates`); the VCF is the `filtered_snps.rsid.vcf` you produced at the end of [the annotation tutorial](../varinat_calling/variant_calling_annotation.md) (Step 6).

> If you do not have your own files, ready-made copies are on the server under  
> `/mnt/nas1/proj/omicss26/ngs_data_analysis/variant_calling/data/`  
> — the BAMs in `bam_clean/`, the filtered VCFs in `vcf/`. The rsIDs are optional here (you will navigate by coordinate); add them with the one `bcftools` command from Step 6 of the annotation tutorial if you want the variants named.

---

## Tools

- `IGV` — the genome browser itself, run on your **own computer**  
    Download it from [igv.org](https://igv.org/doc/desktop/#DownloadPage/); it needs Java 21+ (see the [README](./README.md) for how to check)
- `samtools` — to index a BAM, if its index is missing
- `bgzip` + `tabix` — to compress and index the VCF

> `samtools`, `bgzip` and `tabix` live on the server — run them there, before you copy the files down.

---

## Preparing your files

**Indexes.** IGV cannot open a BAM or a VCF without its index, and the index must sit in the **same folder** as the file it belongs to.

For the BAMs, `MarkDuplicates` already wrote an index when you ran it with `--CREATE_INDEX true` — GATK names it `wes46_dedup.bai`. IGV reads that happily, but the name it recognises most reliably is `wes46_dedup.bam.bai`, which `samtools index` writes:

```bash
samtools index wes46_dedup.bam
```

Note that you index the **BAM**, not the index — this writes `wes46_dedup.bam.bai` next to it.

A VCF needs a different kind of index. Compress it with `bgzip` and index it with `tabix`:

```bash
bgzip -c filtered_snps.rsid.vcf > filtered_snps.rsid.vcf.gz
tabix -p vcf filtered_snps.rsid.vcf.gz
```

Then copy all six files to your computer, open IGV, and set the genome to **hg38** — the same reference you aligned against. Load the two BAMs and the VCF; the two samples should appear one above the other.

> One of these two people is healthy. The other has been diagnosed with **Familial Mediterranean Fever (FMF)**. Nothing in the files tells you which — that is what you are here to work out.

---

## Part 1 — Getting your bearings

### Exercise 1
Type `MEFV` into the navigation box. The view is empty until you do — remember that this is not whole-exome data, so most of the genome has nothing in it, and the reads appear only where there is data.

Once there: **where are the reads?** Describe the shape of the coverage across the gene. Why does it look like this, and not like an even blanket over the whole region?

### Exercise 2
Zoom in on the boundary between a covered region and an empty one. The reads do not stop exactly where the exon does.

Why do reads spill past the ends of exons — and why, in a **whole-exome** experiment, do you nonetheless see variants sitting in introns?

### Exercise 3
Look at the two samples side by side. Do they have the same coverage in the same places? If one sample has reads somewhere the other does not, is that biology, or something else?

---

## Part 2 — Reading the pileup

A read that agrees with the reference is drawn plain grey; a base that disagrees shows up as a coloured tick against it. That is all a SNP is in IGV — colour where the reads and the genome part ways. Zygosity is just the *proportion* of reads that carry the colour.

### Exercise 4
Navigate to **chr16:3,247,073**. Both samples are well covered here, and they disagree.

For each sample: how many reads carry the reference base, and how many carry the alternative? Which sample is **homozygous** and which is **homozygous reference**? How does each of them look in IGV — what is it about the picture that tells you?

### Exercise 5
Now go to **chr16:3,254,573**, where both samples carry the same variant, and it is **heterozygous** in both.

How does a heterozygote look different from the homozygote you just saw? Roughly what fraction of the reads carry the variant, and why that fraction?

### Exercise 6
Both of the above are in *MEFV*, and both are shared between our two people. If a variant appears in the patient **and** in the healthy person, can it be the cause of the disease? What does that let you rule out?

---

## Part 3 — Do you believe the call?

This is the part a VCF cannot answer for you.

### Exercise 7
Go to **chr16:3,243,922**. The caller reports one of the samples as **homozygous for the alternative allele** here — a confident-sounding statement.

Count the reads it is based on. Would you sign your name to that call? What would you want before you did?

> Look at the genotype fields for this variant in your VCF: `1/1:0,2:2:6:49,6,0`. Read off `AD` (how many reads per allele), `DP` (total depth) and `GQ` (how sure the caller is, in Phred). Does `GQ` agree with your own impression?

### Exercise 8
Go to **chr16:3,254,762**. In your VCF, one of the two samples has the genotype `./.` at this position.

Look at that sample in IGV. **What is actually there?** Does `./.` mean the sample carries the reference allele — or something else entirely? Why does this distinction matter when you are trying to prove that someone does *not* have a mutation?

### Exercise 9
Scan along the gene and find a position where you can plainly see a variant in the reads, but the VCF has no entry for it.

Give the coordinates, and explain why the caller passed it over. (There is more than one possible reason — how many can you find? Consider read depth, mapping quality, strand, and the filters you applied in Step 7 of the [pipeline](../varinat_calling/variant_calling_pipeline.md).)

### Exercise 10
Find a variant where the reads supporting it are all in the same orientation, or all crowded at the end of their reads. Why is IGV showing you this, and which of the GATK filters — `FS`, `SOR`, `ReadPosRankSum` — is designed to catch it?

---

## Part 4 — Your candidate

### Exercise 11
Navigate to the FMF candidate you identified in the [annotation tutorial](../varinat_calling/variant_calling_annotation.md).

- How many reads support it, and how many support the reference?
- Do the counts match the `AD` field of the VCF? If not — remember that IGV shows you every read, while the caller only counted the ones that passed its own quality filters. Which reads do you think it threw away?
- Is the variant at the end of the reads, or in the middle? Are both strands represented?

### Exercise 12
Now look at **the other sample** at that same position.

Can you say that this person does not carry the mutation? Be precise about what the evidence supports: is it "no mutation", or "no data"? Look back at Exercise 8.

### Exercise 13
Look at **chr16:3,243,407** — the position of M694V, the most common FMF mutation of all. Neither sample has a variant called here.

What do the reads say? Are both samples homozygous reference at this position, or did nobody have enough coverage to tell? These two possibilities look identical in a VCF and completely different in IGV.

### Exercise 14
Put it together.

**Which of the two samples is the patient, and what is your evidence?**

Then the harder question. FMF is a **recessive** disease — it takes two broken copies of *MEFV* to cause it. Look again at the genotype of the variant you found. Does what you have seen fully explain the diagnosis? If not, what would you go looking for next, and where in this data might it be hiding?

---

## What next

You have now taken the FMF case as far as this dataset allows: from raw reads, through a called and filtered VCF, to an annotated shortlist, and finally back down to the reads to see whether the shortlist tells the truth. Naming the patient — and defending it from the pileup rather than from a line of text — is the whole job a variant analyst does.

Two threads run on from here:

- [**VCF file manipulation**](../vcf_file_manipulation/README.md) — filtering and exploring these same variants at cohort scale with PLINK.
- The **recessive puzzle** from Exercise 14 — a single heterozygous hit does not, on its own, explain a recessive disease. Where might a second variant be hiding, and which of the checks above would you trust to find it?

When everyone has finished, we will compare notes on what you found. Good luck!

#### End of practice ####
