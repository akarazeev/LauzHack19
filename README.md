# LauzHack19 - Service for pathogens and threats identification in hospital environments (by [SOPHiA Genetics](https://www.sophiagenetics.com/home.html))

<img src="img/lauzhack.jpg" href="https://lauzhack.com/" width='20%' align='right'>

As part of the 2019 edition of **LauzHack**, our team has chosen to work on **genomic processing**. This project aims at creating value by substantially streamlining the genomic analysis of patients in a visual and efficient manner.

The way we achieve is by focusing on a pairwise comparison of DNA sequences as well as providing an interactive interface for practitioners.

From setting up a patient's risk profile to visualizating the likelihood of virus infections our end product is powerful, transparent and highly customizable. The clinician will be in control and have an instant diagnosis of his patient.

### Key Features

* Algorithmic
	* <img src="https://biopython.org/assets/images/biopython_logo_s.png" width="100" height="80" />
    * G-Matrix (Patients x Harmful Agents) [@Samuel]
	* Denoising of benign viruses
	* Cross correlation analysis
* Visualization
	* **Tree-of-life** relations
	i.e.: <img src="https://cbsnews1.cbsistatic.com/hub/i/2015/09/21/5816af3f-b563-433f-871c-eb1dabeea491/9efe7bb53aec2efc91071bf1841074b4/tree-of-life-fig1620w.jpg" width="20%"/>
	* Superimposition of histograms
	* Interactive MVP [@Yann21]
	
![](img/heatmap.png)

![](img/hist.png)

### Pipeline
![Flowchart](img/flowchart_backend.png)

```mermaid
graph LR
A(Data selection) -- Data processing --> B{Pairwise Matrix}
B --> C(Additional operations...)
B --> E(Extract parallel information)
E --> D
C --> D(Interface)
```

![Flowchart](img/flowchart_production.png)
```mermaid
graph LR
A(Custom algorithm) --> B(Cloud DB)
B -- Continuous update --> B
B --> D{Hospitals}
D --> E((Clinicians))
C(Front end integration) --> D
```

### Team

* Anton
* Nikita
* Samuel
* Yann

### SOPHiA Partnership

<img src="img/sophia.png" href="https://lauzhack.com/" width='30%' align='right'>

SOPHiA GENETICS combines deep expertise in life sciences and medical disciplines with mathematical capabilities in data computing. Our mission is to bring data analytics solutions to market, to support healthcare professionals by maximizing the power of Data-Driven Medicine. We achieve this mission through the global adoption of SOPHiA artificial intelligence. SOPHiA is built using techniques such as statistical inference, pattern recognition and machine learning. This enables SOPHiA to provide equal benefits to all users, unite experts in a gold standard health tech platform, and motivate expert knowledge sharing for a sustainable impact on future patients.

### References

* Human viruses and associated pathologies: https://viralzone.expasy.org/678
* NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
  - FTP Server: ftp://ftp.ncbi.nlm.nih.gov/refseq/release/
