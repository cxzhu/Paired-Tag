# Paired-Tag protocol

#### To access the latest [Paired-Tag protocol](https://github.com/cxzhu/Paired-Tag/blob/master/protocol/Paired-Tag_protocol.pdf).

### FAQ:
#### 1. What purification method should I use for ordering the oligo DNA sequences?
(a) For 96-well barcode plates, a standard desalting is OK for this purpose.
(b) For Tn5 adaptor sequences (BC Plate#01, DNA_#01_RE to DNA_#12_RE and pMENTs), we used standard desalting. However, some lab reported using HPLC purification for the Tn5 adaptor sequences will increase the tagmentation efficiency; we did not tested this.
(c) for Index PCR primers (N5XX, P7XX and P5 Universal), HPLC purification is required. Note: the "NNNNNNNN" in N5XX and P7XX is the index sequence, not random base mix, for details, please see [Illumina Adaptor Sequences](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html).


#### 2. How can I purify the protein A-Tn5?
The protein A-Tn5 and purification protocol are from Henikoff Lab: 
"pA-Tn5" vector from https://www.addgene.org/124601/. 
The purification is described in [this paper](https://urldefense.com/v3/__https://www.ncbi.nlm.nih.gov/pubmed/31036827__;!!Mih3wA!TjUh5F8Wk_JFXUTsoQYIk0HNossGOompLefbfxxnffxzj8_pxWzxgiXwkMglGLNYPA$).

#### 3. Can I use crosslinked nuclei?
Light crosslinking will help reduce nuclei clumping in [CUT&Tag](https://www.nature.com/articles/s41596-020-0373-x). However, Paired-Tag uses a different library preparation strategy and crosslinking will reduce the DNA library complexity, and thus native nuclei is recommended.

#### 4. How can I reduce nuclei clumping?
Adding 1% BSA will help to reduce nuclei clumping for native nuclei, as described in [this paper](https://www.nature.com/articles/s41587-021-00869-9). According to our recent optimization, adding 1%-5% BSA into Complete Buffer, Med Buffer #1 and Med Buffer #2 and 0.1% BSA into reverse transcriptione buffer will help reduce nuclei clumping.
