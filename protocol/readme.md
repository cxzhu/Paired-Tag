# Paired-Tag protocol & FAQs

#### To access the latest [Paired-Tag protocol](https://github.com/cxzhu/Paired-Tag/raw/master/protocol/Protocol_Github.pdf).

### FAQs:
#### What purification method should I use for ordering the oligo DNA sequences?
(a) For 96-well barcode plates, a standard desalting is OK for this purpose.
(b) For Tn5 adaptor sequences (BC Plate#01, DNA_#01_RE to DNA_#12_RE and pMENTs), we used standard desalting. However, some lab reported using HPLC purification for the Tn5 adaptor sequences will increase the tagmentation efficiency; we did not tested this.
(c) for Index PCR primers (N5XX, P7XX and P5 Universal), HPLC purification is required. Note: the "NNNNNNNN" in N5XX and P7XX is the index sequence, not random base mix, for details, please see [Illumina Adaptor Sequences](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html).

#### What is the purpose of 3'ddC in pMENTs sequence?
The 3'ddC here is to prevent unwanted tailing of pMENTs sequence during TdT tailing process, the tailed oligos may have very long length and cannot be removed from size selection.

#### Whay BC Plate #03 have the names R04-XXX?
The BC Plate #03 is for 2nd round of ligation, the name R04-XXX is to be consistent with our previous Paired-seq protocol. The R03 plate in protocol refers the BC Plate #03 here.

#### Do I need to use secondary antibody?
We don’t use secondary antibody as it may generate fragments with too short sizes, that are not well compatible with the downstream size selection steps. Please double check the primary antibody you plan to use have strongly affinity with protein A. [Protein A wikipedia](https://en.wikipedia.org/wiki/Protein_A#cite_note-11)

#### How can I purify the protein A-Tn5?
The protein A-Tn5 and purification protocol are from Henikoff Lab: 
"pA-Tn5" vector from https://www.addgene.org/124601/. 
The purification is described in [this paper](https://urldefense.com/v3/__https://www.ncbi.nlm.nih.gov/pubmed/31036827__;!!Mih3wA!TjUh5F8Wk_JFXUTsoQYIk0HNossGOompLefbfxxnffxzj8_pxWzxgiXwkMglGLNYPA$).

#### Can I use crosslinked nuclei?
Light crosslinking will help reduce nuclei clumping in [CUT&Tag](https://www.nature.com/articles/s41596-020-0373-x). However, Paired-Tag uses a different library preparation strategy and crosslinking will reduce the DNA library complexity, and thus native nuclei is recommended.

#### What is the rational of the P5 adaptor mix? What are the functions of P5-FokI and P5H-FokI?
The FokI cutting may generate a C in the first base of Read1. To increase sequencing quality, we mix P5H (add one additional H just after Read1 primer) to unshift this consistent C base to increase sequencing quality.

#### How can I reduce nuclei clumping?
Adding 1% BSA will help to reduce nuclei clumping for native nuclei, as described in [this paper](https://www.nature.com/articles/s41587-021-00869-9). According to our recent optimization, adding 1%-5% BSA into Complete Buffer, Med Buffer #1 and Med Buffer #2 and 0.1% BSA into reverse transcriptione buffer will help reduce nuclei clumping.

#### Can I use Nextera XT DNA Library Prep Kit for cDNA library prepration?
Yes, you can use Nextera XT DNA Library Prep Kit to replace Tn5-Adaptor tagmentation. But please note to use TruSeq i7 primer + Nextera i5 primer to do the index PCR amplification, as the cellular barcodes were attached with a TruSeq i7 adaptor during ligation.

#### Are the RNA-seq library stranded?
The RNA-seq library are stranded and the RNA reads transcribed from each strand can be seperated.

#### How to ensure the tagmentation step to prepare RNA library won’t chop away the P7 adapter/cell barcode?
It is hard to ensure tagmentation step do not chop P7/cell barcode. Fortunately, the tagmentation step is carried out in pre-amplified product, so we only need to find a cDNA/Tn5 ratio that majority of the tagmentation events do not cut the adaptors and use size selection to remove the too short fragments.
According to our experience (and our batches of Tn5 protein), use 1 µL 0.01mg/mL Tn5-Adaptor A is good for 10-100 ng of pre-amplified (and NotI digested) products. You may need to do a pilot titration test for your own Tn5.
An alternative way is to use the NextEra XT kit (https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/nextera-xt-dna.html), the Tn5 protein provided in this kit have a wider range of input material. You can use the tagmentation buffer, Tn5 protein, quenching buffer and PCR mix from them, but you will need to use your own P7/N5 primers (not the N7/N5 NextEra primers associated with this product).

#### What are the recommended on sequencing depth?
We typically sequence the libraries to ~25k DNA reads + 25k RNA reads per nuclei. Different histone modifications may have very different library complexities, I would recommend to sequence both DNA and RNA libraries to ~40-60% PCR duplication ratio.
