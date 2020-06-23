library(bookdown)
args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 
cfold=as.character(args[1])#first argument is the number of slave
print(paste0("R will compile in",cfold))
folder="~/public_html/report_lastonly/"
render_book(input="index.Rmd",output_dir="~/public_html/report_lastonly/")
