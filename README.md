# Human CD8 TIL atlas
This repository is describing the construction of a scRNA-seq human T cell ProjecTILs reference. This reference aims at being comprehensive/exhaustive but also interpretable and robust. 

- To construct our model, training data were obtained from filtered data derived from Nick Borcherding collection: https://github.com/ncborcherding/utility

- CD8+ T cells were isolated using scGate. 

- Cycling cells and ISG+ cells were removed using UCell.

- To improve integration, semi-supervised STACAS (ssSTACAS) was run using pre-defined scGate classes as prior knowledge.

The goal for this reference serve as a basis to understand currently overlooked parallels between mice and human settings, and serve as a consistent reference to name clusters across studies.

## Prior knowledge / scGate
<img src="https://user-images.githubusercontent.com/34238952/218471799-633bd162-c581-4324-bb17-099a3ec6dedd.png" height = "400">

     


## Semi-supervised integration

Integration using ssSTACAS https://carmonalab.github.io/STACAS.demo/STACAS.demo.html#semi-supervised-integration
![Capture d’écran 2023-02-03 à 15 21 24](https://user-images.githubusercontent.com/34238952/216626818-af97baa7-d4d9-4de9-a04a-62e8befdb5d3.png)


# Tutorial to use this reference:

Detecting Tpex with this reference:
- https://carmonalab.github.io/ProjecTILs.demo/Detecting_Tpex_CD8_reference.html

Recovering cell types from cell states clusters (eg, activation, cycling, IFN) with this reference:
- https://carmonalab.github.io/ProjecTILs.demo/Transient_gene_programs.html

For more projection examples, you can find them on this repo:
