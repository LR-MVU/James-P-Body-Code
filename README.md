# James-P-Body-Code
This program produces additional analysis on P-Body and mRNA information in the dendrites and soma of neurons. It requires the user to input folder paths for the following
  - Dendritic and Soma prints
  - Soma P-Body results
  - Dendrite P-Body results
  - Soma FISH-QUANT results
  - Dendrite FISH-QUANT results

The user must set these paths within the code itself. An example of these inputs can be found in the sample data folder (https://mcgill.sharepoint.com/:f:/s/VeraUgaldeLab_Group/EsMnrurWy1lLocodaA9DACUByWuXotHNKZm8-tLmC5nfgw?e=Y94hyh). If the user inputs data that does not match across all streams the program prints an error message with information on data read in from each stream (specifically a list containing lists of image names and number of dendrites in image). This allows the user to use a text compare software to delete files that do not match. The sample data folder runs for the dendrites and produces the excel file with results but demonstrates the error throwing process for the soma.

The results spreadsheet produced by this program contains the image name, the dendrite number, the fraction of total dendrite area taken up by P-Bodies, number of P-Bodies per pixel of dendrite area, number of mRNA per pixel of dendrite area, number of mRNA per pixel of P-Body area, and the fraction of total mRNA contained inside the P-Bodies.
