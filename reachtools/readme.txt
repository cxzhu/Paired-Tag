
Update 220524:

1. Merge determine_type into combine2/3 function. Now the readType file should be generated after calling combine2/3, 
so there is no more need to call determine_type separately. This is also to prevent setting wrong barcode length 
when calling determine_type separately.

========================================================================================
Update 220404:

1. Add combine3 to process 8-bp barcodes;
2. Add determine_type and filtBam to mark and remove cross-modality contamination reads;
3. Add duplicates marking function to rmdup2.

========================================================================================

As the library structure of Paired-tag and current version of Paired-seq is different from initial version of Paired-seq described in Zhu et.al., NSMB, 2019, please be aware to use the right version of reachtools.

This one is for Paired-tag/Paired-seq version 2020 only!
For analysis of Paired-seq version 2019 data, please use https://github.com/cxzhu/Paired-seq.

