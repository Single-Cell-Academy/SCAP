# **News**

We have added this News tab to the SCAP app. Come back here to find the newest updates and functions we have implemented into SCAP!

### **2021-07-14 (0.2.2) : Fixed bugs that caused gene names or metadata to be turned into numbers in certain situations. ** 

### **2021-06-24 (v0.2.0 + v0.2.1) : Major update to file conversion module and fixed an issue with annotation.**

### **2021-06-23 (v0.1.2): Updated docker and internal R version used from 4.0.4 to 4.1.0**

### **2021-06-21 (v0.1.1): Fixed minor annotation bug**

### **2021-05-26 : Major update to version 0.1.0**

In this update, we have implemented a number of new features:
* Differential expression calculation between groups
* Enabled loading of .h5ad files for additional modalities alongside RNA
* Internal upgrade from R version 3.6 to 4.0.4

### **2021-03-25 : New data format**

We have completely modified the input format used by SCAP. Instead of using custom loom files, SCAP now uses h5ad formatted files, similar to cellxgene. This means that the data used as input for SCAP can seamlessly be imported into cellxgene as well. We made this decision in an effort for easier compatibility of SCAP analysis with other tools. If you have an existing loom file from previous analysis, you first have to transform this into an h5ad file using the File conversion tab.


### **2020-10-21 : Nebulosa and Scibet added**

We are excited to add two new features to SCAP with this update. 

1) You can now use [Nebulosa](https://bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/introduction.html#:~:text=Nebulosa%20is%20an%20R%20package,convolution%E2%80%9D%20of%20the%20cell%20features) to visualize expression of multiple features on your embeddings. Try it out in the Main tab.

2) We have added [SciBet](http://scibet.cancer-pku.cn/document.html) as a tool to rapidly annotate your data against one of > 90 prepared reference datasets. To find more details about SciBet and try out this new feature, head over to the SciBet tab and give it a try!