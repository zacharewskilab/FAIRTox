#source("./global.R")
source("./heatmap_expr_prep.R")
source("./heatmap_binary_prep.R")

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

hm_pal <- c("#0000FF","#0202FF","#0505FF","#0707FF","#0A0AFE",
           "#0C0CFF","#0F0FFF","#1111FE","#1414FF","#1616FF",
           "#1919FF","#1C1CFF","#1E1EFF","#2121FF","#2323FF",
           "#2626FF","#2828FF","#2B2BFE","#2D2DFF","#3030FF",
           "#3333FF","#3535FF","#3838FF","#3A3AFF","#3D3DFF",
           "#3F3FFF","#4242FF","#4444FF","#4747FF","#4949FE",
           "#4C4CFF","#4F4FFF","#5151FF","#5454FF","#5656FF",
           "#5959FF","#5B5BFF","#5E5EFF","#6060FF","#6363FF",
           "#6666FF","#6868FF","#6B6BFF","#6D6DFF","#7070FF",
           "#7272FF","#7575FF","#7777FF","#7A7AFF","#7C7CFF",
           "#7F7FFF","#8282FF","#8484FF","#8787FF","#8989FF",
           "#8C8CFF","#8E8EFF","#9191FF","#9393FF","#9696FF",
           "#9999FF","#9B9BFF","#9E9EFF","#A0A0FF","#A3A3FF",
           "#A5A5FF","#A8A8FF","#AAAAFF","#ADADFF","#AFAFFF",
           "#B2B2FF","#B5B5FF","#B7B7FF","#BABAFF","#BCBCFF",
           "#BFBFFF","#C1C1FF","#C4C4FF","#C6C6FF","#C9C9FF",
           "#CCCCFF","#CECEFF","#D1D1FF","#D3D3FF","#D6D6FF",
           "#D8D8FF","#DBDBFF","#DDDDFF","#E0E0FF","#E2E2FF",
           "#E5E5FF","#E8E8FF","#EAEAFF","#EDEDFF","#EFEFFF",
           "#F2F2FF","#F4F4FF","#F7F7FF","#F9F9FF","#FCFCFF",
           "#FFFFFF",
           "#FFFCFC","#FFFAFA","#FFF8F8","#FEF6F6","#FFF4F4",
           "#FFF2F2","#FEF0F0","#FFEEEE","#FFECEC","#FFEAEA",
           "#FFE8E8","#FFE6E6","#FFE4E4","#FFE2E2","#FFE0E0",
           "#FFDEDE","#FEDCDC","#FFDADA","#FFD8D8","#FFD6D6",
           "#FFD4D4","#FFD2D2","#FFD0D0","#FFCECE","#FFCCCC",
           "#FFC9C9","#FFC7C7","#FFC5C5","#FEC3C3","#FFC1C1",
           "#FFBFBF","#FFBDBD","#FFBBBB","#FFB9B9","#FFB7B7",
           "#FFB5B5","#FFB3B3","#FFB1B1","#FFAFAF","#FFADAD",
           "#FFABAB","#FFA9A9","#FFA7A7","#FFA5A5","#FFA3A3",
           "#FFA1A1","#FF9F9F","#FF9D9D","#FF9B9B","#FF9999",
           "#FF9696","#FF9494","#FF9292","#FF9090","#FF8E8E",
           "#FF8C8C","#FF8A8A","#FF8888","#FF8686","#FF8484",
           "#FF8282","#FF8080","#FF7E7E","#FF7C7C","#FF7A7A",
           "#FF7878","#FF7676","#FF7474","#FF7272","#FF7070",
           "#FF6E6E","#FF6C6C","#FF6A6A","#FF6868","#FF6666",
           "#FF6363","#FF6161","#FF5F5F","#FF5D5D","#FF5B5B",
           "#FF5959","#FF5757","#FF5555","#FF5353","#FF5151",
           "#FF4F4F","#FF4D4D","#FF4B4B","#FF4949","#FF4747",
           "#FF4545","#FF4343","#FF4141","#FF3F3F","#FF3D3D",
           "#FF3B3B","#FF3939","#FF3737","#FF3535","#FF3333")

DRE <- read.table('./tzheatmap_files/Mouse_DREs_MSS_0.856.txt', sep = '\t', header = TRUE)
AhR <- read.table('./tzheatmap_files/AhR_Binding_FDR_0.05.txt', sep = '\t', header = TRUE)

prepAhR <- function(ahr.dre.input){
  for (g in 1:nrow(ahr.dre.input)){
    if (ahr.dre.input[g,1] %in% DRE[,1]){
      ahr.dre.input[g,2] = 1
      if (ahr.dre.input[g,1] %in% AhR[,1]){
        ahr.dre.input[g,3] = 1
      }
    }
  }
  AHR = heatmap_binary_prep(ahr.dre.input, title = 'ahr-dre', binColumns = c(2,3))
  return(AHR)
}

prepData <- function(gene_list, df, fc_cutoff){
  genes.of.interest <- gene_list
  df.fcCols = c(2:(fc_cutoff+1))
  print(df.fcCols)
  colnames(df[,df.fcCols])
  df.ptCols = c((fc_cutoff+2):ncol(df))
  colnames(df[,df.ptCols])
  DF <- heatmap_expr_prep(df, title = 'dose-response', exprColumns = df.fcCols, pColumns = df.ptCols, fcPal = hm_pal)
  return(DF)
}

heatmap_combine <- function(gene_list, tc, dr = NULL, ahr = NULL, circ = NULL){
  Combined <- rbind(tc, dr, ahr, circ)
  
  #Rescale & relabel Data
  Combined[which(Combined$value < -2 & Combined$type == 'expression'), "value"] = -2
  Combined[which(Combined$value > 2 & Combined$type == 'expression'), "value"] = 2
  Combined[which(Combined$value < 500 & Combined$type == 'count'), "value"] = 500
  Combined[which(Combined$value > 10000 & Combined$type == 'count'), "value"] = 10000
  
  Combined[which(Combined$value == 1 & Combined$type == 'rythm1'), "format"] = "Y"
  Combined[which(Combined$value == 0 & Combined$type == 'rythm1'), "format"] = "N"
  
  Combined[which(Combined$value == 1 & Combined$type == 'rythm2'), "format"] = "X"
  Combined[which(Combined$value == 0 & Combined$type == 'rythm2'), "format"] = ""
  
  Combined[which(Combined$value >= 0.8 & Combined$type == 'significance'), "format"] = "*"
  Combined[which(Combined$value < 0.8 & Combined$type == 'significance'), "format"] = ""
  
  #Reorder Panels and genes
  Combined$Gene = factor(Combined$Gene, levels = rev(gene_list))
  Combined$title = factor(Combined$title, levels = c('ahr-dre','time-course','dose-response','circadian_data'))
  #Combined$variable = factor(Combined$variable, levels = c("pDRE","AhR",colnames(tc.data)[tc.fcCols], colnames(dr.data)[dr.fcCols], 'Counts','rythmicity','rythmChange'))
  
  return(Combined)
}

create_heatmap <- function(Combined){
  hmap <- ggplot(Combined, aes(variable, Gene)) + facet_grid(~ title, scales = 'free_x', space = 'free_x') + theme_minimal() + theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    geom_tile() +
    geom_tile(aes(fill = value), data = Combined %>% filter(type %in% c('expression')), color = 'black', size = 1) +
    scale_fill_gradientn(colours = hm_pal, limits = c(-2,2)) +
    geom_text(aes(label = format), data = Combined %>% filter(type %in% c('significance')), size = 8) +
    
    new_scale('fill') +
    geom_tile(aes(fill = value), data = Combined %>% filter(type %in% c("count")), color = 'black', size = 1) +
    scale_fill_gradient2("Max Count", limits = c(500,10000), low = '#fffe3f', mid = '#fffe3f', high = '#ff3397', midpoint = 3000) +
    
    new_scale('fill') +
    geom_tile(aes(fill = factor(value)), data = Combined %>% filter(type %in% c('binary')), color = 'black', size = 1) +
    scale_fill_manual(values = c("#828282","#00cd42"), name = 'Yes/No', guide = guide_legend(reverse=TRUE)) +
    geom_text(aes(label = format), data = Combined %>% filter(type %in% c('rythm2')), size = 6, color = 'yellow') +
    
    new_scale('fill') +
    geom_tile(fill = 'white', data = Combined %>% filter(type %in% c('rythm1')), color = 'black', size = 1) +
    geom_text(aes(label = format), data = Combined %>% filter(type %in% c('rythm1'))) +
    
    
    scale_x_discrete(position = "top") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 0, color = 'black'),
      axis.text.y = element_text(size = 12, face = "italic", color = 'black'),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  return(hmap)
}
