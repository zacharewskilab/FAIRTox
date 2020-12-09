#' build_dataframes.R
#' Author: Jack Dodson
#' Michigan State University - Zacharewski Lab

builddfDoseResponse <- function(annotation, symbol, connection){
  # Edge cases
  if(symbol == ""){
    return(data.frame(x1 = c()))
  }
  if(length(symbol) > 1){
    stop("More than one gene selected.")
  }
  # Query the database
  sqlcmd <- paste0("SELECT E.Project_ID, Longname, Dose, FoldChange, Species_common, Sex, Organ_name, 
                    P1t, PMID, GEO, DesignType, Chemical_Name, Assay_Name, Strain_Name FROM 
                    (SELECT GeneUID FROM annotation WHERE upper(annotation.'",annotation,"') == upper('",symbol,"')) as A
                    JOIN StudyDesigns as X on X.StudyDesign_ID = E.StudyDesign_ID
                    JOIN ExpressionChange as E ON E.GeneUID = A.GeneUID
                    JOIN FAIR as F ON F.Project_ID = E.Project_ID
                    JOIN Project as P ON P.Project_ID = E.Project_ID
                    JOIN Species as S ON S.Species_ID = E.Species_ID
                    JOIN Sex as Sx ON Sx.Sex_ID = E.Sex_ID
                    JOIN Organ_name as O ON O.Organ_name_ID = E.Organ_name_ID
                    JOIN Chemical as C on C.Chemical_ID = E.Chemical_ID
                    JOIN Assays as Q on Q.Assay_ID = E.Assay_ID
                    JOIN Strain on Strain.Strain_ID = E. Strain_ID
                    ORDER BY E.Project_ID ASC, Dose ASC;")
  df <- data.frame(dbGetQuery(connection, sqlcmd))
  
  # Modify column types
  df$Dose <- as.numeric(as.character(df$Dose))
  df$FoldChange <- as.numeric(as.character(df$FoldChange))
  df$P1t <- as.numeric(as.character(df$P1t))
  df$Project_ID <- as.character(as.numeric(df$Project_ID))
  
  return(df)
}

builddfTimeCourse <- function(annotation, symbol, connection){
  # Edge cases
  if(symbol == ""){
    return(data.frame(x1 = c()))
  }
  if(length(symbol) > 1){
    stop("More than one gene selected.")
  }
  # Query the database
  sqlcmd <- paste0("SELECT E.Project_ID, Longname, TimePoint, FoldChange, Species_common, Sex, Organ_name, 
                    P1t, PMID, GEO, DesignType, Chemical_Name, Assay_Name, Strain_Name FROM 
                   (SELECT GeneUID FROM annotation WHERE upper(annotation.'",annotation,"') == upper('",symbol,"')) as A
                   JOIN StudyDesigns as X on X.StudyDesign_ID = E.StudyDesign_ID
                   JOIN ExpressionChange as E ON E.GeneUID = A.GeneUID
                   JOIN FAIR as F ON F.Project_ID = E.Project_ID
                   JOIN Project as P ON P.Project_ID = E.Project_ID
                   JOIN Species as S ON S.Species_ID = E.Species_ID
                   JOIN Sex as Sx ON Sx.Sex_ID = E.Sex_ID
                   JOIN Organ_name as O ON O.Organ_name_ID = E.Organ_name_ID
                   JOIN Chemical as C on C.Chemical_ID = E.Chemical_ID
                   JOIN Assays as Q on Q.Assay_ID = E.Assay_ID
                   JOIN Strain on Strain.Strain_ID = E. Strain_ID
                   ORDER BY E.Project_ID ASC, Dose ASC;")
  df <- data.frame(dbGetQuery(connection, sqlcmd))
  
  # Modify column types
  df$TimePoint <- as.numeric(as.character(df$TimePoint))
  df$FoldChange <- as.numeric(as.character(df$FoldChange))
  df$P1t <- as.numeric(as.character(df$P1t))
  df$Project_ID <- as.character(as.numeric(df$Project_ID))
  
  return(df)
}

builddfCircadianTop <- function(annotation, symbol, connection){
  # Edge cases
  if(symbol == ""){
    return(data.frame(x1 = c()))
  }
  if(length(symbol) > 1){
    stop("More than one gene selected.")
  }
  # Query the database
  sqlcmd <- paste("SELECT E.Project_ID, ZT, FoldChange, Longname, P1t, Organ_name, DesignType FROM 
                  (SELECT GeneUID FROM annotation WHERE upper(annotation.'",annotation,"') == upper('",symbol,"')) as A
                  JOIN StudyDesigns as X on X.StudyDesign_ID = E.StudyDesign_ID
                  JOIN ExpressionChange as E ON E.GeneUID = A.GeneUID
                  JOIN Project as P ON P.Project_ID = E.Project_ID
                  JOIN Organ_name as O ON O.Organ_name_ID = E.Organ_name_ID
                  ORDER BY E.Project_ID ASC, DOSE ASC;", sep = "")
  df <- data.frame(dbGetQuery(connection, sqlcmd))
  
  # Modify column types
  df$FoldChange <- as.numeric(as.character(df$FoldChange))
  df$P1t <- as.numeric(as.character(df$P1t))
  df$Project_ID <- as.character(as.numeric(df$Project_ID))
  df$ZT <- (as.numeric(df$ZT))
  df$Longname <- as.character(df$Longname)
  df$Organ_name <- as.character(df$Organ_name)
  df$DesignType <- as.character(df$DesignType)
  return(df)
}

builddfCircadianBottom <- function(annotation, symbol, connection){
  # Edge cases
  if(symbol == ""){
    return(data.frame(x1 = c()))
  }
  if(length(symbol) > 1){
    stop("More than one gene selected.")
  }
  # Query the database
  sqlcmd <- paste("SELECT ZT, Project_ID, Dose, NormalizedSignal FROM (SELECT * FROM Animal WHERE Project_ID == 150) as A 
                  JOIN Biosample ON Biosample.Animal_ID = A.Animal_ID 
                  JOIN Animal_treatment ON Animal_treatment.AnimalTreatment_ID = A.AnimalTreatment_ID 
                  JOIN Normalized_Expression ON Normalized_Expression.Biosample_ID = Biosample.Biosample_ID
                  JOIN Species ON Species.Species_ID = A.Species_ID
                  JOIN Sex ON Sex.Sex_ID = A.Sex_ID
                  JOIN Organ ON Organ.Organ_ID = Biosample.Organ_ID
                  JOIN Organ_name ON Organ_name.Organ_name_ID = Organ.Organ_name_ID
                  JOIN Chemical ON Chemical.Chemical_ID = Animal_treatment.Chemical_ID
                  JOIN Assays ON Assays.Assay_ID = Biosample.Assay_ID
                  JOIN Strain ON Strain.Strain_ID = A.Strain_ID
                  JOIN Annotation ON Annotation.GeneUID = Normalized_Expression.GeneUID 
                  WHERE annotation.'",annotation,"' == '",symbol, "';", sep = "")
  df <- data.frame(dbGetQuery(connection, sqlcmd))
  
  # Modify column types
  df$Dose <- as.factor(df$Dose)
  df$ZT <- as.factor(df$ZT)
  df$NormalizedSignal <- 2^df$NormalizedSignal
  return(df)
}

#' Query the database to produce a data frame for plotting functions
#' 
#' @param annotation The annotation type (Symbol, Ensembl ID, or NCBI ID)
#' @param symbol The gene symbol to search for (e.g. Cyp1a1)
#' @param connection The connection to the database
#' @return a data frame produced from querying the dbZach database
#' @examples
#' builddfBE(input$annotationInput, input$SymbolInput, dbconn)
builddfBE <- function(annotation, symbol, connection){
  # Edge cases
  if(symbol == ""){
    return(data.frame(x1 = c()))
  }
  if(length(symbol) > 1){
    stop("More than one gene selected.")
  }
  # Query the database
  sqlcmd <- paste("SELECT * FROM Normalized_Expression as N 
                  JOIN Biosample ON Biosample.Biosample_ID = N.Biosample_ID 
                  JOIN Animal ON Animal.Animal_ID = Biosample.Animal_ID 
                  JOIN Animal_treatment ON Animal_treatment.AnimalTreatment_ID = Animal.AnimalTreatment_ID
                  JOIN Species ON Species.Species_ID = Animal.Species_ID
                  JOIN Sex ON Sex.Sex_ID = Animal.Sex_ID
                  JOIN Organ ON Organ.Organ_ID = Biosample.Organ_ID
                  JOIN Organ_name ON Organ_name.Organ_name_ID = Organ.Organ_name_ID
                  JOIN Chemical ON Chemical.Chemical_ID = Animal_treatment.Chemical_ID
                  JOIN Assays on Assays.Assay_ID = Biosample.Assay_ID
                  JOIN Strain ON Strain.Strain_ID = Animal.Strain_ID 
                  JOIN Annotation ON Annotation.GeneUID = N.GeneUID
                  WHERE Animal_treatment.Dose == 0 AND annotation.'",annotation,"' == '",symbol, "';", sep = "")
  df <- data.frame(dbGetQuery(connection, sqlcmd))
  # Modify column types
  df$Dose <- as.factor(df$Dose)
  df$ZT <- as.factor(df$ZT)
  df$NormalizedSignal <- 2^df$NormalizedSignal
  df$BodyWeight <- as.character(df$BodyWeight)
  df$Lab_Identifier <- as.character(df$Lab_Identifier)
  df$Wet_Weight <- as.character(df$Wet_Weight)
  df$Dry_Weight <- as.character(df$Dry_Weight)
  df$Species_common <- as.character(df$Species_common)
  df$Sex <- as.character(df$Sex)
  df$Strain_Name <- as.character(df$Strain_Name)
  return(df)
}

#' Query the database to produce a data frame for UpSet diagrams
#' 
#' @param selected The projects to analyze
#' @param connection The connection to the database
#' @return a data frame produced from querying the dbZach database
#' @examples
#' buildUpSetdf(input$projectsSelected, dbconn)
buildUpSetdf <- function(selected, connection){
  if(length(selected) < 2){
    stop("Select at least two projects")
  }
  selectedString <- toString(sprintf("'%s'", selected))
  sql_fmt <- paste("SELECT Longname, P1t, Dose, ZT, FoldChange, Symbol, Organ_name FROM
                    (SELECT * from ExpressionChange) as E
                   JOIN Annotation as A on A.GeneUID = E.GeneUID
                   JOIN Organ_name as O on O.Organ_name_ID = E.Organ_name_ID
                   JOIN Project as P on P.Project_ID = E.Project_ID
                   WHERE Longname IN (%s);
                    ", sep = "")
  sqlcmd <- sprintf(sql_fmt, selectedString)
  df <- dbGetQuery(connection, sqlcmd)
  df$Dose <- as.factor(df$Dose)
  df$ZT <- as.factor(df$ZT)
  return(data.frame(df))
}

#TODO: write function signature
buildNewHeatmapdf <- function(selected, connection){
  selectedString <- toString(sprintf("'%s'", selected))
  sql_fmt <- paste0("SELECT Symbol, Dose, TimePoint, Project_ID, FoldChange, P1t, Species_common, Sex, Organ_name FROM (SELECT * FROM ExpressionChange) as E
                    JOIN Annotation as A on A.GeneUID = E.GeneUID
                    JOIN Species as Sp on Sp.Species_ID = E.Species_ID
                    JOIN Sex as Sx on Sx.Sex_ID = E.Sex_ID
                    JOIN Organ_name as O on O.Organ_name_ID = E.Organ_name_ID
                    WHERE Symbol in (%s);")
  sql_cmd <- sprintf(sql_fmt, selectedString)
  df <- dbGetQuery(connection, sql_cmd)
  return(data.frame(df))
}

#' Query the database to produce a dataframe for GSEA
#' 
#' @param selected The genes to analyze
#' @param connection The connection to the database
#' @return A dataframe produced from querying the dbZach database
#' @examples 
#' buildGSEAdf(input$GLD_paste_list_input, dbconn)
buildGSEAdf <- function(selected, connection){
  if(length(selected) == 0){
    stop("Select at least one gene")
  }
  selectedString <- toString(sprintf("'%s'", selected))
  sql_fmt <- paste0("SELECT Symbol, P1t, Dose, FoldChange, Longname, Species_common,
                     Sex, Organ_name, Chemical_Name, Assay_Name, Strain_Name, E.Project_ID FROM
                    (SELECT * from ExpressionChange) as E
                    JOIN Annotation as A on A.GeneUID = E.GeneUID
                    JOIN Species as Sp on Sp.Species_ID = E.Species_ID
                    JOIN Sex as S on S.Sex_ID = E.Sex_ID
                    JOIN Organ_name as O on O.Organ_name_ID = E.Organ_name_ID
                    JOIN Chemical as C on C.Chemical_ID = E.Chemical_ID
                    JOIN Assays as Q on Q.Assay_ID = E.Assay_ID
                    JOIN Strain as St on St.Strain_ID = E.Strain_ID
                    JOIN Project as P on P.Project_ID = E.Project_ID;")
  sql_cmd <- sprintf(sql_fmt, selectedString)
  df <- dbGetQuery(connection, sql_cmd)
  return(data.frame(df))
}

#TODO: write function signature
buildBatchAnalysisdf <- function(connection){
  batchanalysis <- (dbGetQuery(dbconn, "SELECT DISTINCT H.Ensembl_ID, H.H_ID, Project, AID, Species, Strain, Sex, Tissue, Chemical, Dose, Timepoint, ZT, Count  
                                        FROM (SELECT * FROM temp_raw_counts WHERE Count > 0) as R 
                                        JOIN (SELECT * FROM temp_homology) as H ON R.Ensembl_ID = H.Ensembl_ID;"))
  return(batchanalysis)
}
