# Human CD8 TIL atlas
This repository is describing the construction of a scRNA-seq human T cell [ProjecTILs](https://github.com/carmonalab/ProjecTILs) reference. This reference aims at being comprehensive/exhaustive but also interpretable and robust. 

- To construct our model, training data were obtained from filtered data derived from Nick Borcherding [collection](https://github.com/ncborcherding/utility).

- CD8+ T cells were isolated using [scGate](https://github.com/carmonalab/scGate).

- Cycling cells and ISG+ cells were removed using [UCell](https://github.com/carmonalab/UCell).

- To improve integration, [semi-supervised STACAS](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html#semi-supervised-integration) was run using pre-defined scGate classes as prior knowledge.

This ProjecTILs reference map is a ressource that enables rapid cluster annotation across many studies and patients, ultimately decreasing cell annotation time and allowing to draw parallels between studies, diseases and species. 

## Prior knowledge / scGate
To help integrating many datasets, scGate classes for the main CD8 T subsets were designed according to classical immunology markers. The graph below disaplys the classes used, with key marker genes for each cell type.
<br><br/>
<img src="https://user-images.githubusercontent.com/34238952/218471799-633bd162-c581-4324-bb17-099a3ec6dedd.png" height = "400">
<br><br/>
The full list of scGate markers used can be found [here](https://github.com/carmonalab/scGate_models/tree/master/human/CD8_TIL).

## Semi-supervised integration (ssSTACAS)

Integration was performed using [ssSTACAS](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html#semi-supervised-integration) and scGate as prior annotation to help integration. This is especially critical for subsets which are hard to resolve, such as Progenitors Exhausted (TPEX).
![Capture d’écran 2023-02-03 à 15 21 24](https://user-images.githubusercontent.com/34238952/216626818-af97baa7-d4d9-4de9-a04a-62e8befdb5d3.png)


# Tutorial to use this reference:

### Detecting Tpex
- https://carmonalab.github.io/ProjecTILs.demo/Detecting_Tpex_CD8_reference.html

### Recovering cell types from cell states clusters (eg, cycling, IFN, tissue-residency)
- https://carmonalab.github.io/ProjecTILs.demo/Transient_gene_programs.html

### Uncovering activation signal in Yost et al.
- https://carmonalab.github.io/human_CD8_TIL_CaseStudies/Activation_Yost.html

For more projection examples, you can find them here: [ProjecTILs_CaseStudies](https://github.com/carmonalab/ProjecTILs_CaseStudies).
