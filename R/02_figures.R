#' @import ggplot2

#' frequency_Nglycan_components()
#' @description This function will read one or more Nglycan_lists and generate frequency plots for the N-Glycan components
#' @param X can be either a character vector containing Nglycan names OR a list of multiple character vectors, the name of the list components will be use to generate the legend
#' @param type can be either <composition>,<class>,or <both> to indicate the type of ontology terms to be considered in the figure
#' @return a barplot showing the ontology of the selected
#' @author Dusan Velickovic, Geremy Clair
#' @export
freq_ontologies<-function(X,type=c("composition","class","both")){
if(missing(type)){type<-"both"}
  if(is.character(X)){
  l<-list()
  l[[1]]<-X
  X<-l
  names(X)<-"X"
  }

  df<-data.frame(matrix(ncol=4,nrow=0))
  colnames(df)<-c("Ontology","Count","Percentage","Group")
  for (i in 1:length(X)){
    t<-NGlycan_ontologies(X[[i]])
    s<-summary(as.factor(t$Ontology_term))
    s<-data.frame(s)
    s<-cbind(rownames(s),s)
    colnames(s)<-rbind("Ontology","Count")
    s$Percentage<-s$Count/length(X[[i]])*100
    rownames(s)<-NULL
    s$Group <-names(X)[i]
    df<-rbind(df,s)
  }

  unique_o<-unique(df$Ontology)
  unique_g<-unique(df$Group)
  all_combinations <- expand.grid(Ontology = unique_o, Group = unique_g)
  df<-merge(all_combinations,df,by=c("Ontology","Group"),all.x = T)
  df[is.na(df)]<-0

if(type=="class"){
  df<-df[!grepl("Number_of_",df$Ontology),]
  df<-df[order(as.character(df$Ontology)),]
  }
if(type=="composition"){
  df<-df[grepl("Number_of_",df$Ontology),]
  df<-df[order(as.character(df$Ontology)),]
    }
if(type=="both"){
  df1<-df[!grepl("Number_of_",df$Ontology),]
  df1<-df1[order(as.character(df1$Ontology)),]
  df2<-df[grepl("Number_of_",df$Ontology),]
  df2<-df2[order(as.character(df2$Ontology)),]
  df<-rbind(df1,df2)
}

  df$Ontology <- factor(df$Ontology, levels = unique(df$Ontology))

  p<- ggplot(df, aes(x=Ontology,y=Percentage,fill=Group)) +
      geom_bar(stat="identity",position=position_dodge(width = 0.9)) +
      xlab("Ontology Terms") +
      ylab("Percentage of Glycans") +
      ggtitle("Percentage of Nglycans attributed to each ontology term.") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle =90, hjust = 1))
  p
  }
