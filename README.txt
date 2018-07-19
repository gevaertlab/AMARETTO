## **Pancancer Driver Discovery using AMARETTO. ** ##
## **Capturing pancancer genetically and epigenetically driver genes.**##

This repository contains the code accompanying the manuscript "Module analysis captures pancancer genetically and epigenetically deregulated cancer driver genes for smoking and antiviral response".
We have developed an algorithm called AMARETTO to identify pancancer driver genes. AMARETTO integrates pancancer DNA copy number, DNA methylation and gene expression data into modules to identify pancancer driver genes. 

The algorithm: 
	Step 1 identifies candidate cancer driver genes with tumor-specific DNA copy number or DNA methylation alterations compared to normal tissue. We used GISTIC to identify significantly and recurrently deleted or amplified regions in the genome. Similarly, we used MethylMix to identify recurrently hyper-or hypomethylated genes. 
	Step 2 identifies cancer driver genes by modeling the relationship between (epi)genomic and transcriptomic data on an individual gene basis.
	Step 3 uses the cancer driver genes identified from step 2 and now takes a global approach by dissecting global gene expression data into modules of co-expressed genes. Each module also has an associated regulatory program that connects the cancer driver genes from step 1 with their downstream targets. This regulatory program is modeled using linear regression with elastic net regularization.
 
AMARETTO is an R package that can be installed as a package, or the raw source code can be used. AMARETTO supports downloading and processing TCGA data from Firehose (see https://gdac.broadinstitute.org/).