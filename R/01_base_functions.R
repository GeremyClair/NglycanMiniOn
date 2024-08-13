#' Nglycan_db
#' This character vector contains the list of all the Nglycans currently present in Metaspace (as of 07/12/24).
#' @format A character
#'
"Nglycan_db"

#' NGlycan_miner()
#' @description This function will read a list of N_glycan names located in a character vector (each glycan is one element of the char_vector) the algorithm contained in this function will parse the information contained in the glycan names to define the appartenance of the different glycan to glycan ontology terms.
#' @param X a character vector with each element being the name of a given N-glycan.
#' @return This function returns a table providing the appartenance of N glycan to different glycan ontology terms.
#' @author Dusan Velickovic, Geremy Clair
#' @export
NGlycan_miner<-function(X){
  if (!is.character(X)){stop("The imput of the function NGlycan_miner has to be a character vector")}
  t<-data.frame(matrix(nrow = length(X), ncol=27))
  colnames(t)<-c("ID","Hex", "HexNAc", "dHex", "NeuAc","NeuGc","HexA","Su","Ph","Kdn","Me","Ac","Pent","Complex_without_epitope","Polylactosamine","Tetraantenary","Bisecting", "Scialic","Fucosylated","Monofucosylated", "Multifucosylated","Pauci-mannose","Hybrid","High-mannose","Plant","Non-Human","Ambigous_bisecting_or_tetraantenary")
  t$ID<-X

  extract_comp_count<-function(name,comp){
    pattern<-paste0("\\b", comp, ":\\d+\\b")
    match<-regexpr(pattern,name)
    if (match[1] == -1) {return(NA)}
    comp_match<-regmatches(name,match)
    comp_count<-as.integer(sub(paste0(comp,":"),"",comp_match))
    return<-comp_count
  }

  for (i in 2:13){
    t[,i]<-sapply(X,extract_comp_count,comp=colnames(t)[i])
  }

# Now we will establish the appartenance of specific N_glycans to given Ontologies using a set of rules
  #Complex rule: ((HexNAc= (4 or 5) and Hex>2) or (HexNAc=5 and Hex=3)) and (all other=0)
  t$Complex_without_epitope<-(!is.na(t$HexNAc)&!is.na(t$Hex)&(t$HexNAc==4|t$HexNAc==5) &t$Hex>2) &
             (!(t$HexNAc==5 & t$Hex==5)) &
              (is.na(t$dHex) & is.na(t$NeuAc) & is.na(t$NeuGc) & is.na(t$HexA) & is.na(t$Su) & is.na(t$Ph) & is.na(t$Kdn) & is.na(t$Me) & is.na (t$Ac) & is.na(t$Pent))
  #Polylactosamine rule: Hex>7 and HexNAc>6
  t$Polylactosamine <- (!is.na(t$HexNAc)&!is.na(t$Hex)&t$HexNAc>6&t$Hex>7)
  #Tetraantennary rule: HexNAc≥6
  t$Tetraantenary <- (!is.na(t$HexNAc) & t$HexNAc>=6) & (!(t$Hex==6&t$HexNAc==6))
  #Bisecting rule: Hex=HexNAc= (5 or >6)
  t$Bisecting <- (!is.na(t$Hex)&!is.na(t$HexNAc)) & (t$Hex==t$HexNAc) & (t$Hex==5 | t$Hex>6)
  #Sialic rule: NeuAc>0
  t$Scialic <- !is.na(t$NeuAc) & t$NeuAc>0
  #Fucosylated rule: dHex>=1
  t$Fucosylated <- !is.na(t$dHex) & t$dHex>=1
  #Monofucosylated rule: dHex=1
  t$Monofucosylated <- !is.na(t$dHex) & t$dHex==1
  #Polyfucosylated rule: dHex>1
  t$Multifucosylated <- !is.na(t$dHex) & t$dHex>1
  #Pauci-Mannose rule: HexNAc<3
  t$`Pauci-mannose` <- !is.na(t$HexNAc) & !is.na(t$Hex) & t$HexNAc>0 & t$HexNAc<3 &t$Hex<4
  #Hybrid rule:(HexNAc=3 and Hex>4) or (HexNAc>3 and Hex=2)
  t$Hybrid <- (!is.na(t$HexNAc)& !is.na(t$Hex)) & ((t$HexNAc==3 & t$Hex>4)|(t$HexNAc>3 & t$Hex==2))
  #high mannose rule: Hex>3 and HexNAc=(1 or 2) and dHex≥0 and all others =0
  t$`High-mannose`<- (!is.na(t$HexNAc)&!is.na(t$Hex) & t$Hex>3 & (t$HexNAc==1| t$HexNAc==2)) &
    (is.na(t$NeuAc) & is.na(t$NeuGc) & is.na(t$HexA) & is.na(t$Su) & is.na(t$Ph) & is.na(t$Kdn) & is.na(t$Me) & is.na (t$Ac) & is.na(t$Pent))
  #Plant rule: Pent>0
  t$Plant <- !is.na(t$Pent) & t$Pent>0
  #Non-Mammal rule: NeuGc>0
  t$`Non-Human` <- !is.na(t$NeuGc) & t$NeuGc>0
  #Ambigous_bisecting_or_tetraantennary rule: Hex=HexNAc=6
  t$Ambigous_bisecting_or_tetraantenary <- !is.na(t$Hex) & !is.na(t$HexNAc) & t$Hex==t$HexNAc & t$Hex==6
  return(t)
  }

#' NGlycan_ontologies()
#' @description This function will generate ntologies from a list of N-Glycan Names. To get the results in a tabular format instead of a list of ontology terms, please use the function NGlycan_miner
#' @param X a character vector with each element being the name of a given N-glycan.
#' @return This function returns a list of ontology terms associated with the different NGlycans in the list and provide their association with the N-Glycans
#' @author Dusan Velickovic, Geremy Clair
#' @export
NGlycan_ontologies<-function(X){
  if (!is.character(X)){stop("The imput of the function NGlycan_miner has to be a character vector")}
  t<-NGlycan_miner(X)
  o<-data.frame(matrix(ncol=2,nrow=0))
  colnames(o)<-c("Nglycan_name","Ontology_term")

  for(i in 2:13){
    a<-data.frame(cbind(t[,1],t[,i]))
    colnames(a)<-c(colnames(t)[1],colnames(t)[i])
    a<-a[!is.na(a[,2]),]
    if(nrow(a)>0){
    a[,2]<- paste0("Number_of_",colnames(t)[i],"=",a[,2])
    colnames(a)<-c("Nglycan_name","Ontology_term")
    o<-rbind(o,a)
    }
  }

  for(i in 14:ncol(t)){
    a<-data.frame(cbind(t[,1],t[,i]))
    colnames(a)<-c(colnames(t)[1],colnames(t)[i])
    a<-a[a[,2]==TRUE,]
    if(nrow(a)>0){
    a[,2]<-colnames(a)[2]
    colnames(a)<-c("Nglycan_name","Ontology_term")
    o<-rbind(o,a)
    }
  }
o<-o[order(o[,1]),]
return(o)
}

#' Nglycan_Fisher()
#' @description Allows to perform a FisherExact test for N-Glycan lists (query vs universe), by default the universe is the list of all observed N-glycans in Metaspace
#' @param query should be a character vector containing N-Glycans names, the N-glycans in this list have to be present in the universe list
#' @param universe should be a character vector containing N-Glycans names
#' @return a data.frame containing enrichment results.
#' @author Geremy Clair
#' @export
Nglycan_Fisher<-function(query, universe){
  if (!is.character(query)){stop("The 'query' has to be a character vector.")}
  if (missing(universe)){warning("The 'universe' was missing, the Nglycan_db was used as universe.")
    universe<-Nglycan_db}
  if (!is.character(universe)){stop("The 'universe' has to be a character vector.")}
  if (sum(query %in% universe)<length(query)){
    warning("All the elements of the 'query' need to be part of the 'universe', if you are using the default universe, it is possible that one or more of your Nglycans are not currently part of the Nglycan_db.")
    warning("the Nglycans not present in the universe were added to it")
    universe<-c(universe,query[!query %in% universe])
  }

  print(paste0("Your 'query' contains ",length(query)," Nglycans."))
  print(paste0("Your 'universe' contains ",length(universe)," Nglycans."))

  if(sum(duplicated(query)>0)){
    warning(paste0("Your 'query' contained ", sum(duplicated(query)), " duplicated Nglycan name(s)."))
    warning(paste0(sum(duplicated(query)), " duplicated name was removed from the query."))
    print(query[duplicated(query)])
    query<-query[!duplicated(query)]
    print(paste0("Your 'query' now contains ",length(query)," Nglycans."))
  }

  if(sum(duplicated(universe)>0)){
    warning(paste0("Your 'universe' contained ", sum(duplicated(universe)), " duplicated Nglycan name(s)."))
    warning(paste0(sum(duplicated(universe)), " duplicated name was removed from the universe."))
    print(universe[duplicated(universe)])
    universe<-universe[!duplicated(universe)]
    print(paste0("Your 'universe' now contains ",length(universe)," Nglycans."))
  }

  o<-NGlycan_ontologies(universe)
  o_query<-o[o$Nglycan_name %in% query,]
  unique_terms<- unique(o_query[,2])

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:length(unique_terms)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(o_query[,2]==unique_terms[i])
    contingency_list[[i]][1,2]<- sum(o[,2]==unique_terms[i])
    contingency_list[[i]][2,1]<- length(query)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
  }
  names(contingency_list)<-unique_terms

  #create Final table
  fisher_results<-data.frame(matrix(NA, nrow=length(unique_terms), ncol=12))
  colnames(fisher_results)<- c("Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","Test_p","Test_padj","fold_change","Nglycan_in_query","-log10(p)")
  fisher_results$Term_description<-unique_terms
  fisher_results$Pop_query<-length(query)
  fisher_results$Pop_universe<-length(universe)

  #run the test and populate the table
  for (i in 1:length(unique_terms)){
    test<-fisher.test(contingency_list[[i]])
    fisher_results$Test_p [i]<-test$p.value
    fisher_results$Count_query[i]<-contingency_list[[i]][1,1]
    fisher_results$Count_universe[i]<-contingency_list[[i]][1,2]
    fisher_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query)*100
    fisher_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe)*100
    fisher_results$fold_change[i]<- fisher_results$`%_query`[i]/fisher_results$`%_universe`[i]
    fisher_results$Nglycan_in_query[i]<-paste(o_query$Nglycan_name[o_query$Ontology_term==fisher_results$Term_description[i]],collapse =";")
  }

  fisher_results$Test_padj<-p.adjust(p=fisher_results$Test_p,method = "BH",length(fisher_results$Test_p))
  fisher_results$`-log10(p)`<- (-log(fisher_results$Test_p,10))

  fisher_results<-fisher_results[order(fisher_results$Test_p),]
  return(fisher_results)
}

#' Nglycan_EASE()
#' @description Allows to perform a modified FisherExact test (DAVID's EASE test) for N-Glycan lists (query vs universe), by default the universe is the list of all observed N-glycans in Metaspace
#' @param query should be a character vector containing N-Glycans names, the N-glycans in this list have to be present in the universe list
#' @param universe should be a character vector containing N-Glycans names
#' @return a data.frame containing enrichment results.
#' @author Geremy Clair
#' @export
Nglycan_EASE<-function(query, universe){
  if (!is.character(query)){stop("The 'query' has to be a character vector.")}
  if (missing(universe)){warning("The 'universe' was missing, the Nglycan_db was used as universe.")
    universe<-Nglycan_db}
  if (!is.character(universe)){stop("The 'universe' has to be a character vector.")}
  if (sum(query %in% universe)<length(query)){
    warning("All the elements of the 'query' need to be part of the 'universe', if you are using the default universe, it is possible that one or more of your Nglycans are not currently part of the Nglycan_db.")
    warning("the Nglycans not present in the universe were added to it")
    universe<-c(universe,query[!query %in% universe])
  }

  print(paste0("Your 'query' contains ",length(query)," Nglycans."))
  print(paste0("Your 'universe' contains ",length(universe)," Nglycans."))

  if(sum(duplicated(query)>0)){
    warning(paste0("Your 'query' contained ", sum(duplicated(query)), " duplicated Nglycan name(s)."))
    warning(paste0(sum(duplicated(query)), " duplicated name was removed from the query."))
    print(query[duplicated(query)])
    query<-query[!duplicated(query)]
    print(paste0("Your 'query' now contains ",length(query)," Nglycans."))
  }

  if(sum(duplicated(universe)>0)){
    warning(paste0("Your 'universe' contained ", sum(duplicated(universe)), " duplicated Nglycan name(s)."))
    warning(paste0(sum(duplicated(universe)), " duplicated name was removed from the universe."))
    print(universe[duplicated(universe)])
    universe<-universe[!duplicated(universe)]
    print(paste0("Your 'universe' now contains ",length(universe)," Nglycans."))
  }

  o<-NGlycan_ontologies(universe)
  o_query<-o[o$Nglycan_name %in% query,]
  unique_terms<- unique(o_query[,2])

  #create the Contingencies matrix
  contingency_list<-list()
  for (i in 1:length(unique_terms)){
    contingency_list[[i]]<-matrix(ncol=2,nrow=2)
    contingency_list[[i]][1,1]<- sum(o_query[,2]==unique_terms[i])-1
    contingency_list[[i]][1,2]<- sum(o[,2]==unique_terms[i])
    contingency_list[[i]][2,1]<- length(query)-contingency_list[[i]][1,1]
    contingency_list[[i]][2,2]<- length(universe)-contingency_list[[i]][1,2]
  }
  names(contingency_list)<-unique_terms

  #create Final table
  EASE_results<-data.frame(matrix(NA, nrow=length(unique_terms), ncol=12))
  colnames(EASE_results)<- c("Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","Test_p","Test_padj","fold_change","Nglycan_in_query","-log10(p)")
  EASE_results$Term_description<-unique_terms
  EASE_results$Pop_query<-length(query)
  EASE_results$Pop_universe<-length(universe)

  #run the test and populate the table
  for (i in 1:length(unique_terms)){
    test<-fisher.test(contingency_list[[i]])
    EASE_results$Test_p [i]<-test$p.value
    EASE_results$Count_query[i]<-contingency_list[[i]][1,1]
    EASE_results$Count_universe[i]<-contingency_list[[i]][1,2]
    EASE_results$`%_query` [i]<- contingency_list[[i]][1,1]/length(query)*100
    EASE_results$`%_universe`[i]<- contingency_list[[i]][1,2]/length(universe)*100
    EASE_results$fold_change[i]<- EASE_results$`%_query`[i]/EASE_results$`%_universe`[i]
    EASE_results$Nglycan_in_query[i]<-paste(o_query$Nglycan_name[o_query$Ontology_term==EASE_results$Term_description[i]],collapse =";")
  }

  EASE_results$Test_padj<-p.adjust(p=EASE_results$Test_p,method = "BH",length(EASE_results$Test_p))
  EASE_results$`-log10(p)`<- (-log(EASE_results$Test_p,10))

  EASE_results<-EASE_results[order(EASE_results$Test_p),]
  return(EASE_results)
}

#' Nglycan_hypergeometric()
#' @description Allows to perform an hypergeometric test for N-glycan lists (query vs universe), by default the universe is the list of all observed N-glycans in Metaspace
#' @param query should be a character vector containing N-Glycans names, the N-glycans in this list have to be present in the universe list
#' @param universe should be a character vector containing N-Glycans names
#' @return a data.frame containing enrichment results.
#' @author Geremy Clair
#' @export
Nglycan_hypergeometric<-function(query, universe){
  if (!is.character(query)){stop("The 'query' has to be a character vector.")}
  if (missing(universe)){warning("The 'universe' was missing, the Nglycan_db was used as universe.")
    universe<-Nglycan_db}
  if (!is.character(universe)){stop("The 'universe' has to be a character vector.")}
  if (sum(query %in% universe)<length(query)){
    warning("All the elements of the 'query' need to be part of the 'universe', if you are using the default universe, it is possible that one or more of your Nglycans are not currently part of the Nglycan_db.")
    warning("the Nglycans not present in the universe were added to it")
    universe<-c(universe,query[!query %in% universe])
  }

  print(paste0("Your 'query' contains ",length(query)," Nglycans."))
  print(paste0("Your 'universe' contains ",length(universe)," Nglycans."))

  if(sum(duplicated(query)>0)){
    warning(paste0("Your 'query' contained ", sum(duplicated(query)), " duplicated Nglycan name(s)."))
    warning(paste0(sum(duplicated(query)), " duplicated name was removed from the query."))
    print(query[duplicated(query)])
    query<-query[!duplicated(query)]
    print(paste0("Your 'query' now contains ",length(query)," Nglycans."))
  }

  o<-NGlycan_ontologies(universe)
  o_query<-o[o$Nglycan_name %in% query,]
  unique_terms<-unique(o_query$Ontology_term)


#create contingency matrixes
contingency.matrix<- list()
for (i in 1:length(unique_terms))
{
  contingency.matrix[[i]]<- matrix(NA,ncol=2,nrow=2)
  contingency.matrix[[i]][1,1]<- sum(o_query$Ontology_term==unique_terms[i])
  contingency.matrix[[i]][1,2]<- sum(o$Ontology_term==unique_terms[i])
  contingency.matrix[[i]][2,1]<- length(query)-contingency.matrix[[i]][1,1]
  contingency.matrix[[i]][2,2]<- length(universe)-contingency.matrix[[i]][1,2]
  names(contingency.matrix)[i]<-unique_terms[i]
  }


hyper_results<-data.frame(matrix(NA, nrow=length(unique_terms), ncol=12))
colnames(hyper_results)<- c("Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","Test_p","Test_padj","fold_change","Nglycan_in_query","-log10(p)")

hyper_results$Term_description <- unique_terms

hyper_results$Pop_query<- length(query)
hyper_results$Pop_universe<- length(universe)

for (i in 1:length(unique_terms))
{
  test<- min(1-cumsum(dhyper(0:(contingency.matrix[[i]][1,1]-1),contingency.matrix[[i]][1,2],contingency.matrix[[i]][2,2],contingency.matrix[[i]][2,1])))
  hyper_results$Test_p[i]<-test
  hyper_results$Count_query[i]<- contingency.matrix[[i]][1,1]
  hyper_results$Count_universe [i]<- contingency.matrix[[i]][1,2]
  hyper_results$`%_query`[i]<- contingency.matrix[[i]][1,1]/length(query)*100
  hyper_results$`%_universe`[i]<- contingency.matrix[[i]][1,2]/length(universe)*100
  hyper_results$fold_change [i]<- hyper_results$`%_query`[i]/hyper_results$`%_universe`[i]
  hyper_results$Nglycan_in_query[i]<-paste(o_query$Nglycan_name[o_query$Ontology_term==unique_terms[i]],collapse=";")
}
hyper_results$Test_padj<-p.adjust(hyper_results$Test_p,method = "BH",n = nrow(hyper_results))
hyper_results$`-log10(p)`<- (-1*log(hyper_results$Test_p,10))
hyper_results$Test_p[is.na(hyper_results$Test_p)|hyper_results$Test_p==TRUE|hyper_results$Test_p>1]<-1

return(hyper_results)
}

#' Nglycan_binomial()
#' @description Allows to perform a binomial test for N-glycan lists (query vs universe), by default the universe is the list of all observed N-glycans in Metaspace
#' @param query should be a character vector containing N-Glycans names, the N-glycans in this list have to be present in the universe list
#' @param universe should be a character vector containing N-Glycans names
#' @return a data.frame containing enrichment results.
#' @author Geremy Clair
#' @export
Nglycan_binomial<-function(query, universe){
  if (!is.character(query)){stop("The 'query' has to be a character vector.")}
  if (missing(universe)){warning("The 'universe' was missing, the Nglycan_db was used as universe.")
    universe<-Nglycan_db}
  if (!is.character(universe)){stop("The 'universe' has to be a character vector.")}
  if (sum(query %in% universe)<length(query)){
    warning("All the elements of the 'query' need to be part of the 'universe', if you are using the default universe, it is possible that one or more of your Nglycans are not currently part of the Nglycan_db.")
    warning("the Nglycans not present in the universe were added to it")
    universe<-c(universe,query[!query %in% universe])
  }

  print(paste0("Your 'query' contains ",length(query)," Nglycans."))
  print(paste0("Your 'universe' contains ",length(universe)," Nglycans."))

  if(sum(duplicated(query)>0)){
    warning(paste0("Your 'query' contained ", sum(duplicated(query)), " duplicated Nglycan name(s)."))
    warning(paste0(sum(duplicated(query)), " duplicated name was removed from the query."))
    print(query[duplicated(query)])
    query<-query[!duplicated(query)]
    print(paste0("Your 'query' now contains ",length(query)," Nglycans."))
  }

  o<-NGlycan_ontologies(universe)
  o_query<-o[o$Nglycan_name %in% query,]
  unique_terms<-unique(o_query$Ontology_term)

  #create contingency matrixes
  contingency.matrix<- list()
  for (i in 1:length(unique_terms))
  {
    contingency.matrix[[i]]<- matrix(NA,ncol=2,nrow=2)
    contingency.matrix[[i]][1,1]<- sum(o_query$Ontology_term==unique_terms[i])
    contingency.matrix[[i]][1,2]<- sum(o$Ontology_term==unique_terms[i])
    contingency.matrix[[i]][2,1]<- length(query)-contingency.matrix[[i]][1,1]
    contingency.matrix[[i]][2,2]<- length(universe)-contingency.matrix[[i]][1,2]
    names(contingency.matrix)[i]<-unique_terms[i]
  }

  binomial_results<-data.frame(matrix(NA, nrow=length(unique_terms), ncol=12))
  colnames(binomial_results)<- c("Term_description","Count_query","Pop_query","Count_universe","Pop_universe","%_query","%_universe","Test_p","Test_padj","fold_change","Nglycan_in_query","-log10(p)")

  binomial_results$Term_description <- unique_terms

  binomial_results$Pop_query<- length(query)
  binomial_results$Pop_universe<- length(universe)

  for (i in 1:length(unique_terms)){
    test<-binom.test(contingency.matrix[[i]][1,1],length(query),p = contingency.matrix[[i]][1,2]/length(universe))
    binomial_results$Test_p[i]<-test$p.value
    binomial_results$`%_query`[i]<- contingency.matrix[[i]][1,1]/length(query)*100
    binomial_results$`%_universe`[i]<- contingency.matrix[[i]][1,2]/length(universe)*100
    binomial_results$Count_query[i]<- contingency.matrix[[i]][1,1]
    binomial_results$Count_universe [i]<- contingency.matrix[[i]][1,2]
    binomial_results$`%_query`[i]<- contingency.matrix[[i]][1,1]/length(query)*100
    binomial_results$`%_universe`[i]<- contingency.matrix[[i]][1,2]/length(universe)*100
    binomial_results$fold_change [i]<- binomial_results$`%_query`[i]/binomial_results$`%_universe`[i]
    binomial_results$Nglycan_in_query[i]<-paste(o_query$Nglycan_name[o_query$Ontology_term==unique_terms[i]],collapse=";")
  }
  binomial_results$Test_padj<-p.adjust(binomial_results$Test_p,method = "BH",n = nrow(binomial_results))
  binomial_results$`-log10(p)`<- (-1*log(binomial_results$Test_p,10))
  binomial_results$Test_p[is.na(binomial_results$Test_p)|binomial_results$Test_p==TRUE|binomial_results$Test_p>1]<-1

  return(binomial_results)
}

#' Nglycan_KS()
#' @description Allows to perform Nglycan Ontology enrichment using a Kolmogorov-Smirnov test, in this case, the "rank" of the Nglycan will be used to perform the enrichment the rank can be directly provided or derived from a list of pvalues or fold changes.
#' @param rankingTable should be a two column data.frame containing the Nglycan names as first column and the ranking values as second column
#' @param order should either be "ascending" or "descending" to indicate the order of the ranking values to use.
#' @return a data.frame containing enrichment results from the GO retrieve from UniProt.
#' @author Geremy Clair
#' @export
Nglycan_KS<-function(rankingTable, order= "ascending"){
  rankingTable<-data.frame(rankingTable)
  if(ncol(rankingTable)>2){
    warning("Your 'ranking table' contained more than 2 columns, only the first two columns were kept for analysis")
    rankingTable<-rankingTable[,1:2]
    }
  colnames(rankingTable)<-c("IDs","Ranking_values")
  if(missing(order)){order<-"ascending"}
  if(!order %in% c("ascending","descending")){
  stop("The 'order' should be either 'ascending' or 'descending'.")
  }

  rankingTable$IDs<-as.character(rankingTable$IDs)
  rankingTable$Ranking_values<-as.numeric(rankingTable$Ranking_values)

  if(order=="ascending"){
    rankingTable<-rankingTable[order(rankingTable$Ranking_values),]
  }else{
    rankingTable<-rankingTable[order(-rankingTable$Ranking_values),]
    }

  if(sum(duplicated(rankingTable$IDs)>0)){
    warning(paste0("Your 'IDs' contained ", sum(duplicated(rankingTable$IDs)), " duplicated Nglycan name(s)."))
    warning(paste0(duplicated(rankingTable$IDs), " duplicated names were removed from the table."))
    print(rankingTable$IDs[duplicated(rankingTable$IDs)])
    rankingTable<-rankingTable[!duplicated(rankingTable$IDs),]
    print(paste0("Your 'rankingTable' now contains ",nrow(rankingTable)," Nglycans."))
  }

  rankingTable$Rank<-1:nrow(rankingTable)
  o<-NGlycan_ontologies(rankingTable$IDs)
  unique_terms<- o[!duplicated(o[,2]),2]

  RT_by_o<-list()
  IDs_by_o<-list()
  for (i in 1:length(unique_terms)){
    IDs_by_o[[i]]<-o$Nglycan_name[o$Ontology_term==unique_terms[i]]
    RT_by_o[[i]]<-rankingTable$Rank[rankingTable$IDs %in% IDs_by_o[[i]]]
    IDs_by_o[[i]]<-paste(IDs_by_o[[i]],collapse = ";")
  }
  names(RT_by_o)=names(IDs_by_o)=unique_terms

  p<-numeric()
  for (i in 1:length(unique_terms)){
      p[i]<-ks.test(RT_by_o[[i]],rankingTable$Rank,alternative="greater")$p.value
  }
  names(p)<-unique_terms

  final_table<-data.frame(matrix(ncol=5,nrow= length(unique_terms)))
  colnames(final_table) <- c("Ontology_term", "Test_p", "Test_padj","Number_in_list", "IDs_in_list")
  final_table$Ontology_term<-unique_terms
  final_table$Test_p<-p
  final_table$Test_padj<-p.adjust(p,method = "BH")
  final_table$IDs_in_list<- unlist(IDs_by_o)
  final_table$Number_in_list<-nchar(final_table$IDs_in_list)-nchar(gsub(";","",final_table$IDs_in_list))+1

  final_table<-final_table[order(final_table$Test_p),]
  return(final_table)
}
