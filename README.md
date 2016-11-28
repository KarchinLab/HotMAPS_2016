# HotMAPS

## About

Hotspot Missense mutation Areas in Protein Structures (HotMAPS) detects somatic mutation hotspot regions in 3D protein structures. Missense mutations in cancer can often be difficult to interpret their effect, but missense mutations occurring at hotspots tend to implicate a driving role in cancer. Hotspot regions are identified at an individual residue basis, and can be of varying sizes (e.g., 1, 5, or other number of residues). Protein structures often have multiple protein chains, which may arise from the same gene and thus be mutated at the same residue position. Identical chains are accounted for in the model to prevent errors induced by assuming chains are independent. For example, an annotated mutation in a gene forming a homodimer will always contain the same mutation in both chains. PDB biological assemblies are used when available to best reveal biologically meaningful 3D hotspots.

## Links

* [Documentation Wiki](http://github.com/KarchinLab/HotMAPS/wiki/Home)
* [Installation](http://github.com/KarchinLab/HotMAPS/wiki/Installation)
* <a href="http://github.com/KarchinLab/HotMAPS/wiki/Quick-start">Quick start example</a>
* <a href="http://github.com/KarchinLab/HotMAPS/wiki/Tutorial-(Exome-scale)">Tutorial (Exome scale)</a>
* <a href="http://github.com/KarchinLab/HotMAPS/wiki/MySQL-database">Loading MySQL database</a>
* <a href="http://github.com/KarchinLab/HotMAPS/wiki/Tutorial-(Mapping-mutations)">Tutorial (Mapping mutations to protein structure)</a>

## Releases

* HotMAPS-1.1.1 Bug fixes for non-exome scale analysis
* HotMAPS-1.1.0 Map mutations to protein structure and several other changes
* HotMAPS-1.0.0 Initial release

## Citation

If you use HotMAPS in your research, please cite:

* Tokheim C, Bhattacharya R, Niknafs N, Gygax DM, Kim R, Ryan M, Masica DL, Karchin R (2016) Exome-scale discovery of hotspot mutation regions in human cancer using 3D protein structure Cancer Research Published OnlineFirst April 28, 2016; doi:10.1158/0008-5472.CAN-15-3190

## Availability

Releases can be found on github at

* [http://github.com/KarchinLab/HotMAPS/releases](http://github.com/KarchinLab/HotMAPS/releases)

## Platform

HotMAPS works on **linux** operating systems. You will need python installed. For further installation details see the installation page.

## Support

Please contact the package developer Collin Tokheim (ctokheim at jhu dot edu) with suggestions, questions or bug reports.
