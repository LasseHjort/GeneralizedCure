linux <- any(grepl("Linux|linux", Sys.info()))
if(linux){
  project <- "~/GeneralizedCure/"
  setwd(project)
  
  #Figure and table directories
  fig.out <- "Output/Figures"
  tab.out <- "Output/Tables"
}else{
  #Create directory
  project <- "K:/FORSK-Projekt/Projekter/Scientific Projects/110_PhD_Lasse_2015/Projekter/GeneralizedCure/"
  setwd(project)
  
  #Figure, table, and data directories (for thesis and article)
  thesis <- FALSE
  if(thesis){
    fig.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Thesis/papers/paperD/Output/Figures/"
    tab.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Thesis/papers/paperD/Output/Tables/"
  }else{
    fig.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Generalized cure models - BJ/Output/Figures/"
    tab.out <- "C:/Users/sw1y/Dropbox/Apps/ShareLaTeX/Generalized cure models - BJ/Output/Tables/" 
  }
}

#Set generated data directory
data.out <- "GeneratedData"