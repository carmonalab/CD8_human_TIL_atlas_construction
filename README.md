# Human CD8 TIL atlas
This project is aimed at constructing a scRNA-seq human T cell reference. Our aimed reference should be both comprehensive/exhaustive but also interpretable and robust (ie, not having tens of overlapping classes).

- To construct our model, training data were obtained from filtered data derived from Zheng et al. 2021 - Science. https://www.science.org/doi/10.1126/science.abe6474. and Nick Borcherding collection: https://github.com/ncborcherding/utility

- CD8 T cells were isolated using scGate. Cycling cells, HSP+ and ISG+ cells were removed using UCell.

- To improve integration, semi-supervised STACAS (ssSTACAS) was run using pre-defined scGate classes as prior knowledge (cf cartoon below).

- The reference will serve as a basis to understand currently overlooked parallels between mice and human settings, and serve as a consistent reference to name clusters across studies.

![image](https://user-images.githubusercontent.com/34238952/188863492-1ad691b7-af38-49d0-970b-ee15b850432f.png)

Integration using ssSTACAS https://carmonalab.github.io/STACAS.demo/STACAS.demo.html#semi-supervised-integration
![Capture d’écran 2023-02-03 à 15 21 24](https://user-images.githubusercontent.com/34238952/216626818-af97baa7-d4d9-4de9-a04a-62e8befdb5d3.png)

Tutorial to use this reference:

Detecting Tpex with this reference:
- https://carmonalab.github.io/ProjecTILs.demo/Detecting_Tpex_CD8_reference.html

Recovering cell types from cell states clusters (eg, activation, cycling, IFN) with this reference:
- https://carmonalab.github.io/ProjecTILs.demo/Transient_gene_programs.html



