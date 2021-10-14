---
title: "Multiple gene copies"
Objectives:
- Search for another gene in a combination of genomes. 
- Translating command line skills to your specific research. 
---


For bioinformatics study, you know the best about your data and what you would like to achieve with it. The examples shown here can be translated exactly in very few scenarios. However, getting started with command line tools and building pipelines constitute similar steps. Be sure to follow through manuals, input and output formats, and specifically with parameters required for your analysis.

If you have finished the exercises above and interested in additional hands-on experience with similar workflow, we have added several whole genomes from Xanthomonas species and a reference 16S rRNA sequence in the folder ‘/blue/general_workshop/share/16S’. Script ‘slurm_pipeline.sh’ is also added in the same folder. 


~~~
$ cd /blue/general_workshop/<username>

$ cp -r ../share/16S ./

$ cd 16S

$ ls
~~~
{: .language-bash}

~~~
16S_ref.fas
slurm_pipeline.sh
X_albilineans_GPE_PC73.fasta
X_albilineans_Xa-FJ1.fasta
X_axonopodis_Xac29-1.fasta
X_campestris_pv_musacearum_NCPPB_4379.fasta
X_cassavae_CFBP_4642.fasta
X_citri_LMG_9322.fasta
X_citri_UnB-Xtec2D.fasta
X_euvesicatoria_85-10.fasta
X_gardneri_ATCC_19865.fasta
X_hortorum_B07-007.fasta
X_melonis_CFBP_4644.fasta
X_oryzae_NCPPB_4346.fasta
X_oryzae_pv_oryzicola_BLS256.fasta
X_perforans_91-118.fasta
X_perforans_LH3.fasta
X_pisi_CFBP_4643.fasta
X_populi_CFBP_1817.fasta
~~~
{: .output}

We recommend running commands individually from this script. Do you see any kind of variation in the output? How would you resolve this variation?

