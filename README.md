# GDC-RNA-Seq
A simple and preliminary code for processing and analyzing RNA-seq data form GDC

The process relys on awk for parsing all the samples unstranded information. there are two files for this part countM.awk and countM.sh(both should be downloaded and placed int the same folder). 
First thing to do is to download and extract GDC samples
Then open and edit the path in the countM.sh. Set it for the extracted folder then run the script (countM.sh) in your teminal. 

After that the resulting count matrix can be used for further analyzis with R.
