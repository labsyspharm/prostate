#filepath <- "/Users/aliamiryousefi/Desktop/LSP12601.ptrdim.cys"
library(data.table) ##install again
library(googlesheets4)
library(vioplot)



mpath<-"/Users/aliamiryousefi/Desktop/cyf_files2/"
files<- list.files(mpath)

ldata<- list()
for (i in 1:length(files)){
  #names(ldata)[i]<- strsplit(files[i], "\\.p")[[1]][1]
  ldata[[i]]<- fread(cmd=paste("/Users/aliamiryousefi/github/cysift/src/cysift view -R", paste(mpath, files[i], sep="")))
}
for (i in 1:length(files)){names(ldata)[i]<- strsplit(files[i], "\\.p")[[1]][1]}


#dt<- fread(cmd=paste("/Users/aliamiryousefi/github/cysift/src/cysift view -R", filepath)) #restarting the R, the code is dt<- fread(cmd=paste("cysift view -R", filepath))
#the key thing is that the cysift only take the binary input and only the binary input of the cysift format can be read into the cysift


#then save the list data

save(ldata, file = "listdatacysift2_Phen_Ki.RData")
#and load it after that.


load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki.RData")

# This includes all the data in one file. 
# the ops package could be used to decode the binary data into base ten for the gating


## for example the pflag is (xxx)_10 but the gating in 64 bit case is (yyy)_2 for all the 
#fixed-position channels, if a channel is not gated its always zero. 
#install.packages("bitops")
library(bitops)

#the following returns the binary of the base ten value
binary <- function(x) {
  i <- 0
  string <- numeric(64)
  while(x > 0) {
    string[64 - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  first <- match(1, string)
  string[first:64]
}

text<- try(system("/Users/aliamiryousefi/github/cysift/src/cysift view -H ~/Documents/Collaborative_Projects/Jeremiah/Primary_prostate_CyCIF_TME_manuscript_2023/cys_files/LSP12601.ptrdim.cys", intern = TRUE))
MA<- text[grep("@MA", text)]
mar<- gsub("\tIN", "", gsub("@MA\tID:", "", MA))
markers<- mar
for (i in 1:length(mar)){
  markers[i]<- strsplit(mar[i], ":")[[1]][1]
} 

markers<- markers[length(markers):1]

#then we get the adaptive binary length of the vectors

ad_binary <- function(x) {
  i <- 0
  string <- numeric(64)
  while(x > 0) {
    string[64 - i] <- x %% 2
    x <- x %/% 2
    i <- i + 1 
  }
  first <- match(1, string)
  string[(65-length(markers)):64]
}

##now the markers and the ad_binary are parallel to each other!




#sheet<- as.data.table(googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1pXzLJkaX2wkJI_uOdgXLaapg1CJlo459G84FRY1pJYI/edit?pli=1#gid=0", sheet="Sheet1"))

sheet<- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")



high<- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow=="High")
#removing the ones without any TLS
high<- high[-3]

low<- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow=="Low")
low<- low[-c(1,4)]

# for each of ther TLS numbers (different ones) add points for to the main plot
#j is the range for patinets

pdf(file = "a-numbered-TLS-compact2-unicolour.pdf", width=10.5, height=10.5)
par(mfrow = c(4, 4), mar = c(1, 1, 1, 1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
#13000 is the max range of the points
for (j in 17:29){
  f<- ldata[[j]]
  ff<-f
  so<-sample(dim(f)[1], round(dim(f)[1]/50))
  f<- f[so,]
  if (j == 15){f<- ldata[[j]]}
  name<- names(ldata)[j]
  determine_color <- function(nam) {
    if (nam %in% high) {
      return("#7570b3")
    } else if (nam %in% low) {
      return("#1b9e77")
    } else {
      return("grey")
    }
  }
  main_color <- determine_color(name)
  ltls<- length(unique(ff$tls_id))-1
  uni<- unique(f$tls_id)[-1]
  a<- max(ff$x); b<- max(ff$y)
  par(mar = c((10.66-(b/1500)), (10.66-(a/1500)), (10.66-(b/1500)), (10.66-(a/1500))))
  plot(f$x, f$y, pch=19, cex=0.01, cex.axis=1,  # Axis tick labels larger
       cex.lab=1,   # Axis title labels larger
       cex.main=1.4,    # Main title larger 
       col="lightgrey", main = paste(names(ldata)[[j]]), xlab = paste("TLS counts", length(uni), sep = "="), ylab = "", ylim = range(f$y), xlim =range(f$x), col.main = main_color 
       #,panel.first = grid()
  )  
  if (sum(bitAnd(ff$cflag, 8) == 8) > 0){
    ss<- subset(ff, bitAnd(ff$cflag, 8) == 8)
    s1<-sample(dim(ss)[1], round(dim(ss)[1]/50))
    ss<- ss[s1, ]
    points(ss$x, ss$y, col="lightblue", pch=19, cex=0.01)
  }
  for (i in 1:ltls){
    s<- subset(ff, ff$tls_id==uni[i])
    points(s$x, s$y, col="orange", pch=19, cex=0.01)#rainbow(ltls)[i], pch=19, cex=0.01)
    if (dim(s)[1]>0){text(s$x[1], s$y[1], label = uni[i], cex = 0.7)}
  }
  
}
dev.off()


h.ntls<- numeric(length(high))
c<- 1
for (i in high){l<- ldata[[which(names(ldata)== i )]]
h.ntls[c]<- length(unique(l$tls_id))
c<- c+1
}
h.ntls<- h.ntls-1

l.ntls<- numeric(length(low))
c<- 1
for (i in low){l<- ldata[[which(names(ldata)== i )]]
l.ntls[c]<- length(unique(l$tls_id))
c<- c+1
}
l.ntls<- l.ntls-1



h.mtls<- numeric(length(high))
c<- 1
for (i in high){l<- ldata[[which(names(ldata)== i )]]
mean.vector<- numeric(h.ntls[c])
cc<-1
for (j in unique(l$tls_id)[-which(unique(l$tls_id)==0)]){
  mean.vector[cc]<- sum(l$tls_id==j)
  cc<- 1 + cc
}
h.mtls[c]<- mean(mean.vector)
c<- c+1
}

l.mtls<- numeric(length(low))
c<- 1
for (i in low){l<- ldata[[which(names(ldata)== i )]]
mean.vector<- numeric(l.ntls[c])
cc<-1
for (j in unique(l$tls_id)[-which(unique(l$tls_id)==0)]){
  mean.vector[cc]<- sum(l$tls_id==j)
  cc<- 1 + cc
}
l.mtls[c]<- mean(mean.vector)
c<- c+1
}


# Load required libraries
library(ggplot2)
library(dplyr)

#Fig5A
# Generate TLS Counts graph

pdf("TLS-compare.pdf")
set.seed(1236)
gleason <- rep(c("High", "Low"), each = 13)
value <- c(h.ntls, l.ntls)
data <- data.frame(Gleason = gleason, Value = value)

# Create violin plot with data points overlaid and quantile bars
# Conduct t-test
t_test_result <- t.test(Value ~ Gleason, data = data)
p <- ggplot(data, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "TLS Counts", y = "Counts per Samples") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(0, 50) +
  # Add line at y = 45
  geom_segment(x = 1.1, xend = 1.9, y = 45, yend = 45, linetype = "dashed", color = "blue") +
  # Add p-value text
  annotate("text", x = 1.5, y = 45, label = paste("p-value =", 
                                                  formatC(t_test_result$p.value, digits = 3, format = "f")),
           hjust = 0.5, vjust = -1, color = ifelse(t_test_result$p.value < 0.05, "red", "black")) +
  # Specify colors for 'high' and 'low'
  scale_fill_manual(values = c("High" = "#7570b3", "Low" = "#1b9e77"))

# Display the plot
p

#Fig5. B
# Generate Frequency Counts graph
set.seed(1236)
gleason <- rep(c("High", "Low"), each = 13)
value <- c(h.mtls, l.mtls)
data <- data.frame(Gleason = gleason, Value = value)

# Conduct t-test before creating the plot
t_test_result <- t.test(Value ~ Gleason, data = data)

# Create violin plot with data points overlaid and quantile bars
p <- ggplot(data, aes(x = Gleason, y = Value, fill = Gleason)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = "TLS Size", y = "Average number of cells in TLS") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  ylim(0, 9000) +
  # Add line at y = 7500
  geom_segment(x = 1.1, xend = 1.9, y = 7500, yend = 7500, linetype = "dashed", color = "blue") +
  # Add p-value text
  annotate("text", x = 1.5, y = 7500, label = paste("p-value =", 
                                                    formatC(t_test_result$p.value, digits = 3, format = "f")),
           hjust = 0.5, vjust = -1, color = ifelse(t_test_result$p.value < 0.05, "red", "black")) +
  # Specify colors for 'High' and 'Low'
  scale_fill_manual(values = c("High" = "#7570b3", "Low" = "#1b9e77"))

# Display the plot
p
dev.off()
dev.off()


#Fig5D
#Producing the final panel of the Fig 5. 
#markers[-c(27, 28)]
h.mat<- as.data.frame(matrix(0, length(markers), length(high)))
colnames(h.mat)<- high
rownames(h.mat)<- markers
for (j in high){l<- ldata[[which(names(ldata)== j )]]
h.markprop<- numeric(length(markers))
for (i in markers) {
  sum<-0
  for (k in 1:dim(l)[1]){
    sum<- ad_binary(l$pflag[k])[which(markers== i )] + sum
  }
  h.markprop[which(markers == i)] <- sum/dim(l)[1]
}
h.mat[,which(high == j)]<- h.markprop
print(h.mat)
}
#$$$$$

l.mat<- as.data.frame(matrix(0, length(markers), length(low)))
colnames(l.mat)<- low
rownames(l.mat)<- markers
for (j in low){l<- ldata[[which(names(ldata)== j )]]
l.markprop<- numeric(length(markers))
for (i in markers) {
  sum<-0
  for (k in 1:dim(l)[1]){
    sum<- ad_binary(l$pflag[k])[which(markers== i )] + sum
  }
  l.markprop[which(markers == i)] <- sum/dim(l)[1]
}
l.mat[,which(low == j)]<- l.markprop
print(l.mat)
}
#markers<- markers[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]
save(l.mat, file = "l.mat.RData")
save(h.mat, file = "h.mat.RData")

l.mat<- l.mat[-which(colSums(t(l.mat))==0),]
h.mat<- h.mat[-which(colSums(t(h.mat))==0),]
par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
for (i in rownames(l.mat)){
  # Calculate the percentages
  l_percent <- 100 * as.numeric(l.mat[which(rownames(l.mat)== i), ])
  h_percent <- 100 * as.numeric(h.mat[which(rownames(l.mat)== i), ])
  
  # Combine the data
  combined_data <- c(l_percent, h_percent)
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 13), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("orange", "lightblue"))
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups ,              # Data
             method = "jitter", # Random noise
             pch = 19,          # Pch symbols
             col = 1,           # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE)        # Add it over
  
  t_test_result <- t.test(l_percent, h_percent)
  
  # Add t-test result as subtitle
  p_value <- t_test_result$p.value
  if(p_value < 0.1) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }
}

dev.off()




#now doing the analysis for the TLSs only for the producing the markers


#Fig5D with the TLSs
#Producing the final panel of the Fig 5. 
#markers[-c(27, 28)]
h.mat<- as.data.frame(matrix(0, length(markers), length(high)))
colnames(h.mat)<- high
rownames(h.mat)<- markers
for (j in high){l<- ldata[[which(names(ldata)== j )]]
l<-subset(l, l$tls_id!=0)
if (dim(l)[1] > 0){
  h.markprop<- numeric(length(markers))
  for (i in markers) {
    sum<-0
    for (k in 1:dim(l)[1]){
      sum<- ad_binary(l$pflag[k])[which(markers== i )] + sum
    }
    h.markprop[which(markers == i)] <- sum/dim(l)[1]
  }
  h.mat[,which(high == j)]<- h.markprop
  print(h.mat)
}
}

l.mat<- as.data.frame(matrix(0, length(markers), length(low)))
colnames(l.mat)<- low
rownames(l.mat)<- markers
for (j in low){l<- ldata[[which(names(ldata)== j )]]
l<-subset(l, l$tls_id!=0)
if (dim(l)[1] > 0){
  l.markprop<- numeric(length(markers))
  for (i in markers) {
    sum<-0
    for (k in 1:dim(l)[1]){
      sum<- ad_binary(l$pflag[k])[which(markers== i )] + sum
    }
    l.markprop[which(markers == i)] <- sum/dim(l)[1]
  }
  l.mat[,which(low == j)]<- l.markprop
  print(l.mat)
}
}
#markers<- markers[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]
#save(l.mat, file = "l.mat.RData")

l.mat<- l.mat[-which(colSums(t(l.mat))==0),]
h.mat<- h.mat[-which(colSums(t(h.mat))==0),]
par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
for (i in rownames(l.mat)){
  # Calculate the percentages
  l_percent <- 100 * as.numeric(l.mat[which(rownames(l.mat)== i), ])
  h_percent <- 100 * as.numeric(h.mat[which(rownames(l.mat)== i), ])
  
  # Combine the data
  combined_data <- c(l_percent, h_percent)
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 13), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#7570b3", "#1b9e77"), outline=FALSE)
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups ,              # Data
             method = "jitter", # Random noise
             pch = 19,          # Pch symbols
             col = 1,           # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE)        # Add it over
  
  t_test_result <- wilcox.test(l_percent, h_percent, alternative = "two.sided", exact = FALSE, correct = FALSE)
  
  # Add t-test result as subtitle
  p_value <- t_test_result$p.value
  if(p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }
}

dev.off()



#do the above but for all the TLSs as the datapoint

h.mat<- as.data.frame(matrix(0, length(markers), sum(h.ntls)))
#colnames(h.mat)<- high
rownames(h.mat)<- markers
c<-1
aa<-1
for (j in high){l<- ldata[[which(names(ldata)== j )]]
for (t in 1:h.ntls[c]){
  l<-subset(l, l$tls_id == c(unique(l$tls_id)[-which(unique(l$tls_id)==0)])[t])
  if (dim(l)[1] > 0){
    h.markprop<- numeric(length(markers))
    for (i in markers) {
      sum<-0
      for (k in 1:dim(l)[1]){
        sum<- ad_binary(l$pflag[k])[which(markers== i )] + sum
      }
      h.markprop[which(markers == i)] <- sum/dim(l)[1]
    }
    h.mat[,aa]<- h.markprop
  }
  aa<- aa + 1
  print(h.mat)
  l<- ldata[[which(names(ldata)== j )]]
}
c<- c+1
}

l.mat<- as.data.frame(matrix(0, length(markers), sum(l.ntls)))
#colnames(h.mat)<- high
rownames(l.mat)<- markers
c<-1
aa<-1
for (j in low){l<- ldata[[which(names(ldata)== j )]]
for (t in 1:l.ntls[c]){
  l<-subset(l, l$tls_id == c(unique(l$tls_id)[-which(unique(l$tls_id)==0)])[t])
  if (dim(l)[1] > 0){
    l.markprop<- numeric(length(markers))
    for (i in markers) {
      sum<-0
      for (k in 1:dim(l)[1]){
        sum<- ad_binary(l$pflag[k])[which(markers== i )] + sum
      }
      l.markprop[which(markers == i)] <- sum/dim(l)[1]
    }
    l.mat[,aa]<- l.markprop
  }
  aa<- aa + 1
  print(l.mat)
  l<- ldata[[which(names(ldata)== j )]]
}
c<- c+1
}


l.mat<- l.mat[-which(colSums(t(l.mat))==0),]
h.mat<- h.mat[-which(colSums(t(h.mat))==0),]
par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1)) # Set up the plot layout
for (i in rownames(l.mat)) {
  # Calculate the percentages
  l_percent <- 100 * as.numeric(l.mat[which(rownames(l.mat)== i), ])
  h_percent <- 100 * as.numeric(h.mat[which(rownames(l.mat)== i), ])
  
  # Combine the data
  combined_data <- c(l_percent, h_percent)
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 67), rep("H", 190))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#1b9e77", "#7570b3"), outline=FALSE)
  
  # Adjust the y-axis labels to ensure at least two ticks are printed
  axis(2, at = pretty(range(combined_data), n = round(exp(1/median(combined_data)))+10), labels = paste0(pretty(range(combined_data), n = round(exp(1/median(combined_data)))+10), "%"), las=2)
  #print(mean(combined_data))
  # Generate distinctive colors for each group
  colors <- c(rainbow(26)[1:13], rainbow(26)[14:26])
  
  # Replicate colors for each group
  group_colors <- c(rep(colors[1:13], l.ntls), rep(colors[14:26], h.ntls))
  
  # Plot each data point individually with jitteriness and respective color
  for (j in 1:length(combined_data)) {
    jittered_x <- ifelse(groups[j] == "H", 1, 2) + runif(1, -0.2, 0.2)  # Add jitter to x-coordinate
    points(jittered_x, combined_data[j], pch = 20, col = adjustcolor(group_colors[j], alpha.f = 0.5))
  }
  
  # Perform t-test
  t_test_result <- t.test(l_percent, h_percent)
  
  # Add t-test result as subtitle
  p_value <- t_test_result$p.value
  if (p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }
}






#on the tls only 

#bitand for the value I want for the exlusive detection of the cells, (8 is for the tumors of the andreanes ROIs) so for getting those we need the rows with the bitAnd(l$cflag[i], 8)==8 conditions (i is the row)
#checking for them not being zero for the or if there are multiple conditions
#colouring 
#fdr or bonfronei correction
#building the VP tree for all tumor and find the nearest neighbour to those in C++ 
##getting the help of the codes of Jeremiah as build_table, cell_table, ...
#get the jerry workflow for the analysis


library(FNN)
library(data.table) 
library(googlesheets4)
library(vioplot)
library(bitops)

load("/Users/aliamiryousefi/listdatacysift2.RData")

text<- try(system("/Users/aliamiryousefi/github/cysift/src/cysift view -H ~/Documents/Collaborative_Projects/Jeremiah/Primary_prostate_CyCIF_TME_manuscript_2023/cys_files/LSP12601.ptrdim.cys", intern = TRUE))
MA<- text[grep("@MA", text)]
mar<- gsub("\tIN", "", gsub("@MA\tID:", "", MA))
markers<- mar
for (i in 1:length(mar)){
  markers[i]<- strsplit(mar[i], ":")[[1]][1]
} 

markers<- markers[length(markers):1]

sheet<- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")

high<- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow=="High")
#removing the ones without any TLS
high<- high[-3]

low<- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow=="Low")
low<- low[-c(1,4)]

#marker equabalvazne tof the pflat
mpflag<- numeric(length(markers))
for (i in 0:(length(markers)-1)){
  mpflag[i+1]<- 2^i
}
mpflag<- mpflag[length(mpflag):1]

mpflag<- mpflag[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]
markers<- markers[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]


mat<- as.data.frame(matrix(0, length(markers), length(names(ldata))))
rownames(mat)<- markers
colnames(mat)<- names(ldata)


#
mat<- list()
for (i in c(1:15, 17:21, 24:29)){
  element_name <- paste(names(ldata)[i])
  Tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8)
  Strom<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) != 8 & ldata[[i]]$tls_id != 0)
  if (dim(Strom)[1]!=0){
    nStromaTLSs<- length(unique(Strom$tls_id))
    mat[[element_name]] <- as.data.frame(matrix(0, length(markers), nStromaTLSs))
    rownames(mat[[element_name]])<- markers
    colnames(mat[[element_name]])<- as.character(unique(Strom$tls_id))
    for (k in 1:nStromaTLSs){
      for (j in 1:length(markers)){
        strom <- subset(Strom, bitAnd(Strom$pflag, mpflag[j]) == mpflag[j] & Strom$tls_id == unique(Strom$tls_id)[k])
        data_points <- cbind(Tumor$x, Tumor$y)  # Combine x and y into a matrix
        query_point <- matrix(c(strom$x, strom$y), ncol=2)
        knn_result <- get.knnx(data_points, query_point, k=1)
        nearest_neighbor_index = knn_result$nn.index # Extract the index of the nearest neighbor
        nearest_neighbor <- data_points[nearest_neighbor_index, ]  # Get the coordinates of the nearest neighbor
        diffs <- nearest_neighbor - cbind(strom$x, strom$y)  # Calculate the differences in x and y coordinates between corresponding points
        distances <- sqrt(rowSums(diffs^2)) # Calculate the Euclidean distance for each pair
        mat[[element_name]][j, k] <- mean(distances) # Print the pairwise distances
      }
      print(mat)
    }
    #print(mat)
  }
}


for (i in 1:length(names(mat))){
  for (j in 1:length(colnames(mat[[i]]))){
    colnames(mat[[i]])[j]<- paste(names(mat)[i], colnames(mat[[i]])[j], sep="@")
  }
}

save(mat, file = "mat.RData")

load("mat.RData")
matr<- as.data.frame(mat[[1]])
for (i in names(mat)[-1]){
  matr<- cbind(matr, mat[[i]])
}

mean_distances<- colMeans(matr, na.rm = TRUE)*0.325

ld<- log(mean_distances)
ld[which(ld== min(ld))]<- 0


# Extract names from the mean_distances vector
mean_distances_names <- names(mean_distances)

# Filter names that contain any of the high identifiers
filtered_names <- mean_distances_names[sapply(mean_distances_names, function(name) any(sapply(low, function(h) grepl(h, name))))]

# Subset the mean_distances vector to get only the values whose names match the filtered names
filtered_values <- mean_distances[filtered_names]

# Calculate the mean of these filtered values
median_value.low <- median(filtered_values)

# Extract names from the mean_distances vector
mean_distances_names <- names(mean_distances)

# Filter names that contain any of the high identifiers
filtered_names <- mean_distances_names[sapply(mean_distances_names, function(name) any(sapply(high, function(h) grepl(h, name))))]

# Subset the mean_distances vector to get only the values whose names match the filtered names
filtered_values <- mean_distances[filtered_names]

# Calculate the mean of these filtered values
median_value.high <- median(filtered_values)





#mat<- mat[, -c(16, 22, 23)]
# TLS proximity chart
pdf("a-TLS proximity Chart.pdf")
# Assuming 'ica_breakdown' is your list
high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

# Load necessary library for additional colors
library(RColorBrewer)

# Setup the plotting area
par(mar = c(2, 2, 2, 2))
plot(0, type = "n", xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", axes = FALSE, asp = 1)
title(main="TLS Proximity Chart", font.main=2, col.main="navy")

# Define the data - 26 equal sections
values <- rep(1, 26) # Equal values for equal sections

# Generate random numbers for each section to display inside the pie chart
random_numbers <- sample(1:100, 26, replace = TRUE)

# Frequencies to display outside the pie chart
frequencies <- sample(100:1000, 26, replace = TRUE)

# Define 26 colors
colors <- brewer.pal(n = min(26, 8), name = "Dark2")
colors <- colorRampPalette(colors)(26)

colors[1:length(high)] <- "#7570b3"  # Color for high
colors[(length(high) + 1):(length(high) + length(low))] <- "#1b9e77"  # Color for low
#colors_transparent <- adjustcolor(colors, alpha.f = 0.5)
#colors[(length(high) + length(low) + 1):26] <- "grey"  # Remaining colors as grey

# Calculate the starting and ending angles for each slice
angles_deg <- cumsum(c(0, values / sum(values) * 360))
angles_rad <- angles_deg * pi / 180

# Draw each pie slice as a polygon
for (i in 1:length(values)) {
  x <- c(0, 1.4*cos(angles_rad[i]), 1.4*cos(angles_rad[i + 1]), 0)
  y <- c(0, 1.4*sin(angles_rad[i]), 1.4*sin(angles_rad[i + 1]), 0)
  polygon(x, y, col ="white", border = TRUE, lty = 3, lwd = 0.6)
}




radii <- seq(1, 1.4, by = 0.05) # Define radii for the circles
for (r in radii) {
  symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1)
  if (r == 1.2){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 0.7)}
  if (r == 1.2){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1.6)}
  if (r == 1.3){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 0.7)}
  if (r == 1.4){ symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 2.2)}
}
radiiii <- seq(1.39, 1.401, by = 0.001) # Define radii for the circles
for (r in radiiii) {
  symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1.5)
}
#add ther median values for each group
r.high <- 1+ log(median_value.high)/20
r.low  <- 1+ log(median_value.low)/20
symbols(0, 0, circles = r.high, inches = FALSE, add = TRUE, fg = "#7570b3", lwd = 1.8, lty= 5)
symbols(0, 0, circles = r.low, inches = FALSE, add = TRUE, fg = "#1b9e77", lwd = 2, lty= 3)

# Loop through each section to add the random numbers and frequencies


##soreting the data
tumor.portion<- numeric(length(values))
names(tumor.portion)<- c(high,low)
for (i in 1:13){
  tumor.portion[i]<- sum(bitAnd(ldata[[which(names(ldata)==c(high,low)[i])]]$cflag , 8) == 8)/dim(ldata[[which(names(ldata)==c(high,low)[i])]])[1]
}
hi.portion<- names(tumor.portion)[order(tumor.portion[1:13])]
for (i in 14:26){
  tumor.portion[i]<- sum(bitAnd(ldata[[which(names(ldata)==c(high,low)[i])]]$cflag , 8) == 8)/dim(ldata[[which(names(ldata)==c(high,low)[i])]])[1]
}
lo.portion<- names(tumor.portion[14:26])[order(tumor.portion[14:26])]
high<- hi.portion
low<- lo.portion

for (i in 1:length(values)) {
  x <- c(0, cos(angles_rad[i]), cos(angles_rad[i + 1]), 0)
  y <- c(0, sin(angles_rad[i]), sin(angles_rad[i + 1]), 0)
  polygon(x, y, col = adjustcolor(colors[i], alpha.f = 0.05 + sqrt(sum(bitAnd(ldata[[which(names(ldata)==c(high,low)[i])]]$cflag , 8) == 8)/dim(ldata[[which(names(ldata)==c(high,low)[i])]])[1])), border = TRUE)
}
radiii <- seq(0.99, 1.001, by = 0.001) # Define radii for the circles
for (r in radiii) {
  symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "lightgray", lwd = 1.5)
}
for (i in 1:length(values)) {
  # Convert angles to radians for text positioning
  angle_rad_inner <- (angles_deg[i] + angles_deg[i+1]) / 2 * pi / 180
  
  # Calculate sqrt of the sum of squares of sin and cos values
  
  # Coordinates for the random numbers inside the pie
  x_pos_inner1 <- 0.24 * cos(angle_rad_inner)
  y_pos_inner1 <- 0.24 * sin(angle_rad_inner)
  
  x_pos_inner2 <- 0.65 * cos(angle_rad_inner)
  y_pos_inner2 <- 0.65 * sin(angle_rad_inner)
  
  x_pos_inner3 <- 0.9 * cos(angle_rad_inner)
  y_pos_inner3 <- 0.9 * sin(angle_rad_inner)
  
  # Coordinates for the frequencies just outside the pie
  x_pos_outer <- 1.5 * cos(angle_rad_inner)
  y_pos_outer <- 1.5 * sin(angle_rad_inner)
  
  x_pos_outer1 <- 1.1 * cos(angle_rad_inner)
  y_pos_outer1 <- 1.1 * sin(angle_rad_inner)
  
  x_pos_outer2 <- 1.2 * cos(angle_rad_inner)
  y_pos_outer2 <- 1.2 * sin(angle_rad_inner)
  
  x_pos_outer3 <- 1.3 * cos(angle_rad_inner)
  y_pos_outer3 <- 1.3 * sin(angle_rad_inner)
  
  x_pos_outer4 <- 1.4 * cos(angle_rad_inner)
  y_pos_outer4 <- 1.4 * sin(angle_rad_inner)
  # Add the random numbers inside the pie
  #text(x_pos_outer, y_pos_outer, labels = random_numbers[i], cex = 0.7)
  jj<-i
  i<- c(which(names(ldata)%in%c(high, low)[jj]))
  
  text(x_pos_inner2, y_pos_inner2, labels = length(unique(ldata[[i]]$tls_id))-1, cex = (0.65 + (length(unique(ldata[[i]]$tls_id))-1)/33), font=3 )
  
  # Add the frequencies outside the pie; adjust the alignment based on the angle
  if (angle_rad_inner > pi/2 && angle_rad_inner < 3*pi/2) {
    text(x_pos_inner1, y_pos_inner1, labels = names(ldata)[i], cex = 0.7, font= 2, srt = 180 * angle_rad_inner / pi - 180, adj = 1)
  } else {
    text(x_pos_inner1, y_pos_inner1, labels = names(ldata)[i], cex = 0.7, font= 2, srt = 180 * angle_rad_inner / pi, adj = 0)
  }
  if (names(ldata)[i] %in% names(mat)) {
    text(x_pos_outer, y_pos_outer, labels = dim(mat[[which(names(mat)==names(ldata)[i])]])[2], font= 2, cex = (0.65 + (dim(mat[[which(names(mat)==names(ldata)[i])]])[2])/33))
    text(x_pos_inner3, y_pos_inner3, labels = (length(unique(ldata[[i]]$tls_id))-1) - (dim(mat[[which(names(mat)==names(ldata)[i])]])[2]), cex = (0.65 + ((length(unique(ldata[[i]]$tls_id))-1) - (dim(mat[[which(names(mat)==names(ldata)[i])]])[2]))/33))
    err<- runif(1, -0.1, 0.1)
    for (j in 1:length(as.numeric(ld[which(unlist(strsplit(names(ld), "@"))[c(TRUE, FALSE)] == names(ldata)[i])])/20)){
      err<- runif(1, -0.08, 0.08)
      x_pos_data <- (1+ (as.numeric(ld[which(unlist(strsplit(names(ld), "@"))[c(TRUE, FALSE)] == names(ldata)[i])])/20))[j] * cos(angle_rad_inner+err)
      y_pos_data <- (1+ (as.numeric(ld[which(unlist(strsplit(names(ld), "@"))[c(TRUE, FALSE)] == names(ldata)[i])])/20))[j] * sin(angle_rad_inner+err)
      opacity<- 1 - sqrt((x_pos_data)^2 + (y_pos_data)^2) +1
      points(x_pos_data, y_pos_data, pch=20, col=adjustcolor(colors[jj], alpha.f= opacity^2.5))}
  }
  else{
    text(x_pos_outer, y_pos_outer, 0, cex = 0.65, font= 2)
    text(x_pos_inner3, y_pos_inner3, labels = length(unique(ldata[[i]]$tls_id))-1, cex = (0.65 + (length(unique(ldata[[i]]$tls_id))-1)/33))
  }
  if (i == 23) {
    text(x_pos_outer1, y_pos_outer1, paste(7, "um", " "), cex = 0.7, col="darkgrey", font= 2)
    text(x_pos_outer2, y_pos_outer2, paste(50, "um", " "), cex = 0.75, col="darkgrey", font= 2)
    text(x_pos_outer3, y_pos_outer3, paste(400, "um", " "), cex = 0.75, col="darkgrey", font= 2)
    text(x_pos_outer4, y_pos_outer4, paste(3, "mm", " "), cex = 0.8, col="darkgrey", font= 2)
  }
  i<-jj
}

# Drawing 7 outer circles with adjusted radii to fit the new plot dimensions

dev.off()
dev.off()





# vioplot::vioplot(t(mat[order(rowSums(mat)),]*0.325),las=2,main="mtcars$mpg",col="deepskyblue",notch=TRUE)
# par(mfrow = c(3, 7), mar = c(2.2, 3.4, 2.1, 2.2)) # Set up the plot layout
# 
# # Assuming 'mat' has been scaled by 0.325 already
# max_val <- max(mat) * 0.325
# 
# for (i in rownames(mat[order(rowSums(mat)),])) {
#   lp <- as.numeric(mat[which(rownames(mat) == i), ])
#   
#   boxplot(lp * 0.325, range = 3,
#           main = paste(i, "+", sep = ""), 
#           ylab = "", xlab = "", yaxt = "n", 
#           col = c("lightblue"), outline = TRUE, ylim = c(0, 1500))
#   
#   # Adding jittered data points
#   jittered_x <- rep(1, length(lp)) + runif(length(lp), -0.25, 0.25)
#   points(jittered_x, lp * 0.325, pch = 20, col = adjustcolor("brown", alpha.f = 0.5))
#   points(sort(sum(lp))*0.325/27, pch=15, col="black")
#   
#   # Generate tick marks and labels for y-axis with mu symbol and 100 unit intervals
#   ticks <- seq(0, 1500, by = 250) # Adjust the maximum value as needed based on your data
#   labels <- paste0(ticks, "μ")
#   
#   # Draw the custom y-axis
#   axis(2, at = ticks, labels = labels, las = 2)
# }



##### unsupervised clustering of the TLSs

TLSmarkers<- c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")


for (i in 1:length(names(ldata))){
  LSPid<- as.character(rep(names(ldata)[i], dim(ldata[[i]])[1]))
  ldata[[i]] <- cbind(ldata[[i]], LSPid)
}

patients_with_tls<- names(ldata)[-c(2,5,12)]

t<- data.frame()
for (i in patients_with_tls){
  tlsp<- subset(ldata[[i]], ldata[[i]]$tls_id != 0)
  tlsp<- tlsp[, c(1:36, 135)]
  #lsptlsid<- as.character(rep(i, dim(tlsp)[1]))
  #tlspp<- cbind(tlsp, lsptlsid)
  t<- rbind(t, tlsp)
  print(length(unique(tlsp$tls_id)))
}

tlsP<- character()
for (i in 1:dim(t)[1]){
  tlsP[i]<- paste(t$tls_id[i], t$LSPid[i], sep= "@")
}
t<- cbind(t, tlsP)



#t is the table of the all the cells in the TLSs, make sure you have the mpflag and markers vecotrs and proceed with the heatmap

# 
# ##the same way we gonna make the tc as the TLS complement of the t matrix with the inclusion and exxlusion from the tumor area
# tc<- data.frame()
# for (i in patients_with_tls){
#   tlsp<- subset(ldata[[i]], ldata[[i]]$tls_id == 0)
#   tlsp<- tlsp[, c(1:36, 129)]
#   #lsptlsid<- as.character(rep(i, dim(tlsp)[1]))
#   #tlspp<- cbind(tlsp, lsptlsid)
#   tc<- rbind(tc, tlsp)
#   print(length(unique(tlsp$tls_id)))
# }
# 
# 
# ##getting the decomposition of the tissue without tlss (tc) for Tumor and Stroma area!
# tcT<- subset(tc, bitAnd(tc$cflag, 8) == 8)
# 
# tcS<- subset(tc, bitAnd(tc$cflag, 8) != 8)


text<- try(system("/Users/aliamiryousefi/github/cysift/src/cysift view -H ~/Documents/Collaborative_Projects/Jeremiah/Primary_prostate_CyCIF_TME_manuscript_2023/cys_files/LSP12601.ptrdim.cys", intern = TRUE))
MA<- text[grep("@MA", text)]
mar<- gsub("\tIN", "", gsub("@MA\tID:", "", MA))
markers<- mar
for (i in 1:length(mar)){
  markers[i]<- strsplit(mar[i], ":")[[1]][1]
} 

markers<- markers[length(markers):1]

mpflag<- numeric(length(markers))
for (i in 0:(length(markers)-1)){
  mpflag[i+1]<- 2^i
}
mpflag<- mpflag[length(mpflag):1]

mpflag<- mpflag[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]
markers<- markers[-which(markers %in% c("Hoechst", "Ki67", "CD24", "CD19", "TCF1", "GranzymeB", "CD11b"))]


mat<- as.data.frame(matrix(0, length(markers), length(names(ldata))))
rownames(mat)<- markers
colnames(mat)<- names(ldata)

matT<- mat
matS<- mat

##doing three separate ones for all off tls portion of the positive cells, on and off tumor as well

for (i in 1:29){
  pot<- subset(ldata[[i]], ldata[[i]]$tls_id == 0) #the tumor and nonTLS
  potT<- subset(ldata[[i]], ldata[[i]]$tls_id == 0 & bitAnd(ldata[[i]]$cflag, 8) == 8) #the tumor and nonTLS
  potS<- subset(ldata[[i]], ldata[[i]]$tls_id == 0 & bitAnd(ldata[[i]]$cflag, 8) != 8) #the tumor and nonTLS
  for (j in 1:length(markers)){
    mat[j,i] <- dim(subset(pot, bitAnd(pot$pflag, mpflag[j]) == mpflag[j]))[1] / dim(pot)[1]
    matT[j,i] <- dim(subset(potT, bitAnd(potT$pflag, mpflag[j]) == mpflag[j]))[1] / dim(potT)[1]
    matS[j,i] <- dim(subset(potS, bitAnd(potS$pflag, mpflag[j]) == mpflag[j]))[1] / dim(potS)[1]
  }
  print(mat)
}

# mat<- mat[, -c(16, 22, 23)]
# matT<- matT[, -c(16, 22, 23)]
# matS<- matS[, -c(16, 22, 23)]

###making the boxplots
# for mat

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

mat<- mat[, c(which(colnames(mat) %in% low),which(colnames(mat) %in% high))]
matS<- matS[, c(which(colnames(matS) %in% low),which(colnames(matS) %in% high))]
matT<- matT[, c(which(colnames(matT) %in% low),which(colnames(matT) %in% high))]
# Increase the outer margin to allow space for the main title
par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1), oma = c(0, 0, 3, 0))

# Begin your plotting loop
for (i in rownames(mat)){
  # Calculate the percentages
  percent <- 100 * as.numeric(mat[which(rownames(mat)== i), ])
  
  # Combine the data
  combined_data <- percent
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 13), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#7570b3", "#1b9e77"), outline=FALSE)
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups,
             method = "jitter",
             pch = 19,
             col = 1,
             vertical = TRUE,
             add = TRUE)
  
  t_test_result <- wilcox.test(percent[1:13], percent[14:26], alternative = "two.sided", exact = FALSE, correct = FALSE)
  p_value <- t_test_result$p.value
  
  # Add t-test result as subtitle
  if(p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red", outer = FALSE)
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, outer = FALSE)
  }
}

# Add the main title using mtext, placed in the outer top margin
mtext("Percentages on Tissue off TLS", side = 3, outer = TRUE, line = 1, cex = 1.5)
# Load necessary libraries
library(dplyr)

# Assuming matT is your data frame
# Remove columns that are entirely NaN
matT_clean <- matT %>% select_if(~any(!is.na(.)))

# Adjust the plotting parameters as per the original script
par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1), oma = c(0, 0, 3, 0))

# Loop through the rownames of the cleaned matrix
for (i in rownames(matT_clean)) {
  # Calculate the percentages
  percent <- 100 * as.numeric(matT_clean[which(rownames(matT_clean) == i), ])
  
  # Filter out NaN values from percent if necessary
  percent <- na.omit(percent)
  
  # Combine the data
  combined_data <- percent
  
  # Update groups based on the actual length of combined_data
  midpoint <- length(combined_data) / 2
  groups <- c(rep("L", 13), rep("H", 10))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  if(length(combined_data) > 1) {  # Check if there's more than one data point
    boxplot(combined_data ~ groups, range = 3,
            main = paste(i, "+", sep = ""), 
            ylab = "", xlab = "", yaxt = "n", 
            col = c("#7570b3", "#1b9e77"), outline=FALSE)
    axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las = 2)
    stripchart(combined_data ~ groups,
               method = "jitter",
               pch = 19,
               col = 1,
               vertical = TRUE,
               add = TRUE)
    
    # Perform t-test if possible (when there are enough data points for each group)
    if(length(percent[1:midpoint]) > 0 && length(percent[(midpoint + 1):length(percent)]) > 0) {
      t_test_result <- wilcox.test(percent[1:midpoint], percent[(midpoint + 1):length(percent)],
                                   alternative = "two.sided", exact = FALSE, correct = FALSE)
      p_value <- t_test_result$p.value
      
      # Add t-test result as subtitle
      if(p_value < 0.05) {
        subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
        mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red", outer = FALSE)
      } else {
        subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
        mtext(subtitle, side = 3, line = 0.2, cex = 0.8, outer = FALSE)
      }
    }
  }
}

# Add the main title using mtext, placed in the outer top margin
mtext("Percentages on Tumor off TLS", side = 3, outer = TRUE, line = 1, cex = 1.5)


par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1), oma = c(0, 0, 3, 0))
for (i in rownames(matS)){
  # Calculate the percentages
  percent <- 100 * as.numeric(matS[which(rownames(matS)== i), ])
  
  # Combine the data
  combined_data <- percent
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 13), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#7570b3", "#1b9e77"), outline=FALSE)
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups,
             method = "jitter",
             pch = 19,
             col = 1,
             vertical = TRUE,
             add = TRUE)
  
  t_test_result <- wilcox.test(percent[1:13], percent[14:26], alternative = "two.sided", exact = FALSE, correct = FALSE)
  p_value <- t_test_result$p.value
  
  # Add t-test result as subtitle
  if(p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red", outer = FALSE)
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, outer = FALSE)
  }
}

# Add the main title using mtext, placed in the outer top margin
mtext("Percentages on Stroma off TLS", side = 3, outer = TRUE, line = 1, cex = 1.5)


dev.off()



# 
# ##heatmap
# heatmap.m<- matrix(0, 9, 257)   #prebuilding the matrix
# 
# immune.markers<- c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")
# for (i in 1:length(unique(t$tlsP))){
#   each.tls.total<- subset(t, t$tlsP == unique(t$tlsP)[i])
#   for (j in 1:length(immune.markers)){
#     heatmap.m[j,i] <- dim(subset(each.tls.total, bitAnd(each.tls.total$pflag, mpflag[which(markers==immune.markers[j])]) == mpflag[which(markers==immune.markers[j])]))[1] / dim(each.tls.total)[1]
#   }
# }
# 
# h.dframe<- as.data.frame(heatmap.m)
# colnames(h.dframe) <- unique(t$tlsP)
# rownames(h.dframe) <- immune.markers
# 
# heatmap(h.dframe)
# 
# library(ggplot2)
# library(reshape2)
# 
# # Assume df is your data frame
# # df <- your_dataframe_here
# 
# # Melting the data frame
# df_melted <- melt(heatmap.m)
# 
# # Creating the heatmap
# ggplot(data = df_melted, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "blue", high = "red") +
#   theme_minimal() +
#   labs(x = "Immune Markers", y = "TLSs", fill = "Cells + %")
# 
# a<- 100
# # Creating a DataFrame with 257 rows and 8 columns filled with NA
# df_na <- data.frame(replicate(12, rep(NA, length(unique(t$tlsP)))))
# df_na100 <- data.frame(replicate(12, rep(NA, length(unique(t$tlsP))*a)))
# colnames(df_na) <- c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206", "tlsP")
# colnames(df_na100) <- c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206", "tlsP", "Gleason", "GleaTLSno")
# for (i in 1:length(unique(t$tlsP))){
#   rs<- which(t$tlsP==unique(t$tlsP)[i])
#   df_na[i, c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")]<- apply(t[rs, c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")], 2, median, na.rm = TRUE)
#   df_na[i, 10]<- unique(t$tlsP)[i]
#   df_na[i, 10]<- strsplit(df_na[i,10], "@")[[1]][2]
# }
# 
# 
# cut.off <- 10
# highfre<- as.character()
# highinf<- as.character()
# lowfre<- as.character()
# lowinf<- as.character()
# c<-1
# for (i in high){
#   if(length(unique(ldata[[i]]$tls_id))-1 > cut.off) {highfre[c]<- i; c<- c+1}
# }
# c<-1
# for (i in high){
#   if(length(unique(ldata[[i]]$tls_id))-1 <= cut.off) {highinf[c]<- i; c<- c+1}
# }
# c<-1
# for (i in low){
#   if(length(unique(ldata[[i]]$tls_id))-1 > cut.off) {lowfre[c]<- i; c<- c+1}
# }
# c<-1
# for (i in low){
#   if(length(unique(ldata[[i]]$tls_id))-1 <= cut.off) {lowinf[c]<- i; c<- c+1}
# }
# 
# 
# 
# ##we need to have the df_na, high, highfre, highinf, vectors as the LSP ids of the patients with gleason high, gleason high and high number of TLSs, and gleason high with small number of the TLSs (cut-off 10 TLSs!)
# set.seed(123)
# for (i in 1:length(unique(t$tlsP))){
#   rs<- sample(which(t$tlsP==unique(t$tlsP)[i]), a)
#   df_na100[c((i*a-a+1):(i*a)), c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")] <- t[rs, c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")]
#   df_na100[c((i*a-a+1):(i*a)), 10]<- rep(unique(t$tlsP)[i], a)
#   df_na100[c((i*a-a+1):(i*a)), 10]<- rep(df_na[i,10], a)
#   if(df_na100[(i*a), 10] %in% high) {df_na100[c((i*a-a+1):(i*a)), 11]<- "High"}
#   if(df_na100[(i*a), 10] %in% low) {df_na100[c((i*a-a+1):(i*a)), 11]<- "Low"}
#   if(df_na100[(i*a), 10] %in% highfre) {df_na100[c((i*a-a+1):(i*a)), 12]<- "H-Fre"}
#   if(df_na100[(i*a), 10] %in% highinf) {df_na100[c((i*a-a+1):(i*a)), 12]<- "H-Inf"} 
#   if(df_na100[(i*a), 10] %in% lowfre) {df_na100[c((i*a-a+1):(i*a)), 12]<- "L-Fre"}
#   if(df_na100[(i*a), 10] %in% lowinf) {df_na100[c((i*a-a+1):(i*a)), 12]<- "L-Inf"}
#   #else{df_na100[c((i*a-a+1):(i*a)), 12]<- "missing"}
# }
# 
# #df_na1100 is the random 100 value TLS selected for the markers and patients
# data_for_umap <- df_na100[, c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")]
# 
# #data_for_umap <- df_na[, c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")]
# 
# 
# 
# 
# 
##calculating the spatial diversity

####for having it on multiple immune markers
library(fastICA)
library(ggplot2)
library(gridExtra)  # for arranging ggplots

plots <- list()  # Initialize an empty list to store the plots
ica_total <- 0

for (i in seq(1, 16, 1)) {
  subt <- subset(t, t$tlsP == unique(t$tlsP)[i])
  df <- data.frame(subt)

  # Assign categories based on bitwise conditions
  df$category <- ifelse(bitAnd(df$pflag, mpflag[which(markers == immune.markers[1])]) == mpflag[which(markers == immune.markers[1])], "CD20",
                        ifelse(bitAnd(df$pflag, mpflag[which(markers == "CD4")]) == mpflag[which(markers == "CD4")], "CD4",
                               ifelse(bitAnd(df$pflag, mpflag[which(markers == immune.markers[8])]) == mpflag[which(markers == immune.markers[8])], "CD11c",
                                   ifelse(bitAnd(df$pflag, mpflag[which(markers == immune.markers[5])]) == mpflag[which(markers == immune.markers[5])], "CD8a", "Others"))))

  df$mean_x <- mean(df$x)
  df$mean_y <- mean(df$y)

  ica_result <- fastICA(as.matrix(df[, c("x", "y")]), n.comp = 2)
  ica_trace <- sum(apply(ica_result$S, 2, var))  # Corrected assignment
  ica_total <- ica_total + round(10000 * (ica_trace - 2), 0)

  df$ica_end_x1 <- df$mean_x + ica_result$A[1, 1]
  df$ica_end_y1 <- df$mean_y + ica_result$A[2, 1]
  df$ica_end_x2 <- df$mean_x + ica_result$A[1, 2]
  df$ica_end_y2 <- df$mean_y + ica_result$A[2, 2]

  p <- ggplot(df, aes(x, y)) +
    geom_point(aes(color = category), size = 0.65, alpha = 0.8) +  # Added transparency to points
    geom_segment(aes(x = mean_x, y = mean_y, xend = ica_end_x1, yend = ica_end_y1),
                 arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "darkgreen", linewidth = 0.5) +
    geom_segment(aes(x = mean_x, y = mean_y, xend = ica_end_x2, yend = ica_end_y2),
                 arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "darkgreen", linewidth = 0.5) +
    scale_color_manual(values = c("CD20" = "red", "CD4" = "yellow", "CD11c" = "blue", "CD8a" = "green", "Others" = "grey")) +
    ggtitle(paste(paste("ICA Trace:", round(10000 * (ica_trace - 2), 0)), unique(t$tlsP)[i], sep = " / ")) +
    theme_minimal()

  plots[[length(plots) + 1]] <- p
}

gridExtra::grid.arrange(grobs = plots, ncol = 4, nrow = 4)


###calculating the average spatial ICA for each patients


#exluding the gigantic TLSs of LSP12629 image
t.reserve<- t
ica_all<- numeric(length(unique(t.reserve$LSPid)[-12]))
ica_breakdown<- list()
ica_S<- list()
ica_A<- list()
for (k in unique(t.reserve$LSPid)[-12]){
  t<- subset(t.reserve, t.reserve$LSPid == k)
  ica_total<-0
  ica.vector<- numeric(length(unique(t$tlsP)))
  ica.A.vector<- numeric(length(unique(t$tlsP)))
  ica.S.vector<- numeric(length(unique(t$tlsP)))
  for (i in seq(1, length(unique(t$tlsP)), 1)) {
    subt <- subset(t, t$tlsP == unique(t$tlsP)[i])
    df <- data.frame(subt)
    df$category <- ifelse(bitAnd(df$pflag, mpflag[which(markers==immune.markers[1])]) == mpflag[which(markers==immune.markers[1])], "CD20", "Rest")

    df$mean_x <- mean(df$x)
    df$mean_y <- mean(df$y)

    ica_result <- fastICA(as.matrix(df[, c("x", "y")]), n.comp = 2)
    a<- as.numeric(dist(range(ica_result$X[,1])))/2
    b<- as.numeric(dist(range(ica_result$X[,2])))/2
    c<- (a+b)/2
    ica.S.vector[i]<- pi*a*b
    ica.A.vector[i]<- 4/3*pi*a*b*c
    ica_trace <- sum(apply(ica_result$S, 2, var))
    ica.vector[i]<-round(10000 * (ica_trace - 2), 0)
    ica_total<- ica_total+round(10000 * (ica_trace - 2), 0)
  }
  ica_all[which(unique(t.reserve$LSPid)[-12]==k)] <- ica_total/length(unique(t$tlsP))
  ica_breakdown[[which(unique(t.reserve$LSPid)[-12]==k)]] <- ica.vector
  ica_A[[which(unique(t.reserve$LSPid)[-12]==k)]] <- ica.A.vector
  ica_S[[which(unique(t.reserve$LSPid)[-12]==k)]] <- ica.S.vector
}
names(ica_all)<- unique(t.reserve$LSPid)[-12]
names(ica_breakdown)<- unique(t.reserve$LSPid)[-12]
names(ica_A)<- unique(t.reserve$LSPid)[-12]
names(ica_S)<- unique(t.reserve$LSPid)[-12]
t<- t.reserve
###plotting the results
# Assuming 'ica_breakdown' is your list
high <- c("LSP12601", "LSP12625", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12607", "LSP12611", "LSP12613", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

# Sort each subclass by their first quantile values
high_sorted <- high[order(sapply(ica_breakdown[high], function(x) quantile(x, 0.25)))]
low_sorted <- low[order(sapply(ica_breakdown[low], function(x) quantile(x, 0.25)))]

# Combine and order categories
all_categories_sorted <- c(high_sorted, low_sorted)

# Reorder the list based on the ordered categories
ordered_breakdown <- ica_breakdown[all_categories_sorted]

# Assign specific colors: '#7570b3' for high and '#1b9e77' for low
colors <- ifelse(names(ordered_breakdown) %in% high_sorted, "#7570b3", "#1b9e77")

# Create the boxplot with ordered list, assigned colors, paler background, and no outliers
boxplot(ordered_breakdown, las = 2, col = sapply(colors, adjustcolor, alpha.f = 0.5),
        main = "TLSs Spatial ICA by Gleason grade, Ordered by 1st Quantile",
        xlab = "", ylab = "Trace of Spatial ICA",
        staplewex = 0.5, staplecol = colors, border = colors,
        whisklty = 1, whiskcol = colors, plot = TRUE, outline = FALSE)  # outline=FALSE removes outliers

# Color gradient from red to white to blue
color_palette <- colorRampPalette(c("red", "white", "blue"))(71)  # Generates a gradient with 71 steps for 0 to 70

# Add jittered points to the boxplots with fixed gradient
for(i in seq_along(ordered_breakdown)) {
  points(jitter(rep(i, length(ordered_breakdown[[i]]))), ordered_breakdown[[i]],
         col = color_palette[pmax(pmin(ordered_breakdown[[i]], 70), 0) + 1], pch = 16, cex = 1.0)  # Increased point size
}

# Add a legend for clarity
legend("topright", legend = c("High", "Low"), fill = adjustcolor(c("#7570b3", "#1b9e77"), alpha.f = 0.5), title = "Class")



####ordered

# Combine high and low into one array
all_categories <- c(high, low)

# Calculate 10th percentiles for each category and order them
tenth_percentiles <- sapply(ica_breakdown[all_categories], function(x) quantile(x, 0.1))

# Sort the entire array of categories based on their 10th percentile values
all_categories_sorted <- all_categories[order(tenth_percentiles)]

# Reorder the list based on the overall ordered categories
ordered_breakdown <- ica_breakdown[all_categories_sorted]

# Assign specific colors: '#7570b3' for high and '#1b9e77' for low
colors <- ifelse(names(ordered_breakdown) %in% high, "#7570b3", "#1b9e77")

# Create the boxplot with ordered list, assigned colors, paler background, and no outliers
boxplot(ordered_breakdown, las = 2, col = sapply(colors, adjustcolor, alpha.f = 0.5),
        main = "TLSs Spatial ICA, Ordered by 10th Percentile",
        xlab = "", ylab = "Trace of Spatial ICA",
        staplewex = 0.5, staplecol = colors, border = colors,
        whisklty = 1, whiskcol = colors, plot = TRUE, outline = FALSE)  # outline=FALSE removes outliers

# Color gradient from red to white to blue
color_palette <- colorRampPalette(c("red", "white", "blue"))(71)  # Generates a gradient with 71 steps for 0 to 70

# Add jittered points to the boxplots with fixed gradient
for(i in seq_along(ordered_breakdown)) {
  points(jitter(rep(i, length(ordered_breakdown[[i]]))), ordered_breakdown[[i]],
         col = color_palette[pmax(pmin(ordered_breakdown[[i]], 70), 0) + 1], pch = 16, cex = 1.0)  # Increased point size
}

# Add a legend for clarity
legend("topright", legend = c("High", "Low"), fill = adjustcolor(c("#7570b3", "#1b9e77"), alpha.f = 0.5), title = "Class")


######integrating and doing the test for all


# Assuming 'ica_breakdown' is your list
high <- c("LSP12601", "LSP12625", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12607", "LSP12611", "LSP12613", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")

# Aggregate data for high and low categories
high_data <- unlist(ica_breakdown[high])
low_data <- unlist(ica_breakdown[low])

# Create a new list for the aggregated boxplots
aggregated_data <- list(High = high_data, Low = low_data)

# Colors for the boxplots
box_colors <- c("#7570b3", "#1b9e77")  # Colors for High and Low
pale_colors <- sapply(box_colors, adjustcolor, alpha.f = 0.5)

# Color gradient from red to white to blue for points
color_palette <- colorRampPalette(c("red", "white", "blue"))(71)  # Generates a gradient with 71 steps for 0 to 70

# Perform a one-sided t-test
t_test_result <- t.test(low_data, high_data, alternative = "greater")

# Formatting the p-value text
p_value_text <- sprintf("p-value: %.3f%s", t_test_result$p.value, ifelse(t_test_result$p.value < 0.05, "*", ""))

# Plotting
boxplot(aggregated_data, names = c("High", "Low"), col = pale_colors,
        main = "TLSs Spatial ICA by Gleason Grade, Aggregated by Patients",
        xlab = "Gleason", ylab = "Trace of Spatial ICA",
        staplewex = 0.5, staplecol = box_colors, border = box_colors,
        whisklty = 1, whiskcol = box_colors, outline = FALSE)  # outline=FALSE removes outliers

# Increasing jitter amount
jitter_amount = 0.15  # Increased for more spread

# Adding jittered points with color gradient based on value
points(jitter(rep(1, length(high_data)), amount = jitter_amount), high_data,
       col = color_palette[pmax(pmin(high_data, 70), 0) + 1], pch = 16, cex = 1.2)
points(jitter(rep(2, length(low_data)), amount = jitter_amount), low_data,
       col = color_palette[pmax(pmin(low_data, 70), 0) + 1], pch = 16, cex = 1.2)

# Add a legend for class colors
legend("topright", legend = c("High", "Low"), fill = adjustcolor(c("#7570b3", "#1b9e77"), alpha.f = 0.5), title = "Class")

# Display the p-value on the plot
text(x = 1.5, y = max(c(high_data, low_data)) + 1, labels = p_value_text,
     cex = 1.2, col = "red", font = 2)  # font = 2 for bold





###hImmune overal check for the Gleason high and low

##smapling 100000 datapoint from heach slide to form the csv

setwd("/Users/aliamiryousefi/Desktop/GU-csv2/")
set.seed(1987)
Hset<- data.frame()
for (i in 1:length(high)){
  read<- read.csv(file = paste(high[i], ".csv", sep = "_cell"))
  dim<- dim(read)
  imageid<- rep(high[i], dim[1])
  read <- cbind(read, imageid)
  x<- read$X_centroid
  read$X_centroid <- (x - min(x)) / (max(x) - min(x)) + (i - 1) %/% 4 + 1
  y<- read$Y_centroid
  read$Y_centroid <- (y - min(y)) / (max(y) - min(y)) + (i - 1) %% 4 + 1
  Hset<- rbind(Hset, as.data.frame(read[sample(dim[1], 100000, replace = FALSE), ]))
}
write.csv(Hset, file = "Ghigh.csv", row.names = FALSE)
write.csv(Hset[,-(dim(2)+1)], file = "Ghigh-noID.csv", row.names = FALSE)

Lset<- data.frame()
for (i in 1:length(low)){
  read<- read.csv(file = paste(low[i], ".csv", sep = "_cell"))
  dim<- dim(read)
  imageid<- rep(low[i], dim[1])
  read <- cbind(read, imageid)
  x<- read$X_centroid
  read$X_centroid <- (x - min(x)) / (max(x) - min(x)) + (i - 1) %/% 4 + 1
  y<- read$Y_centroid
  read$Y_centroid <- (y - min(y)) / (max(y) - min(y)) + (i - 1) %% 4 + 1
  Lset<- rbind(Lset, as.data.frame(read[sample(dim[1], 100000, replace = FALSE), ]))
}
write.csv(Lset, file = "Glow.csv", row.names = FALSE)
write.csv(Lset[,-(dim(2)+1)], file = "Glow-noID.csv", row.names = FALSE)



###phenotyping based on the real gates!!!



library(FNN)
library(data.table) 
library(googlesheets4)
library(vioplot)
library(bitops)

load("/Users/aliamiryousefi/listdatacysift2.RData")

text<- try(system("/Users/aliamiryousefi/github/cysift/src/cysift view -H ~/Documents/Collaborative_Projects/Jeremiah/Primary_prostate_CyCIF_TME_manuscript_2023/cys_files/LSP12601.ptrdim.cys", intern = TRUE))
MA<- text[grep("@MA", text)]
mar<- gsub("\tIN", "", gsub("@MA\tID:", "", MA))
markers<- mar
for (i in 1:length(mar)){
  markers[i]<- strsplit(mar[i], ":")[[1]][1]
} 

markers<- markers[length(markers):1]

sheet<- read.csv("/Users/aliamiryousefi/Desktop/GU-material/Prostate SPORE - Sheet1.csv")

high<- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow=="High")
#removing the ones without any TLS
high<- high[-3]

low<- subset(sheet$CyCIF.Slide.ID, sheet$Gleason.HighLow=="Low")
low<- low[-c(1,4)]

#marker equabalvazne tof the pflat
mpflag<- numeric(length(markers))
for (i in 0:(length(markers)-1)){
  mpflag[i+1]<- 2^i
}
mpflag<- mpflag[length(mpflag):1]

immune.markers<- c("CD20", "CD68", "CD163", "CD3d", "CD8a", "FOXP3", "PD1", "CD11c", "CD206")
## now there are two type of phenotyping, one is to get the phenotype list for pehnotypes and their positive and negative, then form a vector of "others", and in that for all the phenotypes,
#start from the most abondant one and fill in overlaying on each other, so the smaller ones will alwasy get the priority of visibility,
#this is to avoid the overblowing of the vector above the lenght of the data. 

## the second method is onlly to get the side bar plots for every thing that is positive showing its frequcny which would be resulting to more than one due to inevitable defficieny with 
#gating strategy to leave the false positives for any phenotypes

##method one/ based on the view of the html of the cyft tool 

#forming the frequncy vector
library(bitops)
###coarse phenotype table
for (i in 1:length(names(ldata))){
  d<- ldata[[i]]
  coarse_phen_vec<- rep("Others", length(d$pflag))
  coarse_phen_vec[bitAnd(d$pflag, mpflag[which(markers=="SMA")]) == mpflag[which(markers=="SMA")] | bitAnd(d$pflag, mpflag[which(markers=="CD44")]) == mpflag[which(markers=="CD44")]] <- "Stromal cells"
  coarse_phen_vec[bitAnd(d$pflag, mpflag[which(markers=="CD31")]) == mpflag[which(markers=="CD31")]] <- "Endothelial cells"
  coarse_phen_vec[bitAnd(d$pflag, mpflag[which(markers=="CD20")]) == mpflag[which(markers=="CD20")]] <- "B cells"
  coarse_phen_vec[bitAnd(d$pflag, mpflag[which(markers=="CD3d")]) == mpflag[which(markers=="CD3d")]] <- "T cells"
  coarse_phen_vec[bitAnd(d$pflag, mpflag[which(markers=="CD11b")]) == mpflag[which(markers=="CD11b")] | bitAnd(d$pflag, mpflag[which(markers=="CD68")]) == mpflag[which(markers=="CD68")] | bitAnd(d$pflag, mpflag[which(markers=="CD11c")]) == mpflag[which(markers=="CD11c")] | bitAnd(d$pflag, mpflag[which(markers=="CD163")]) == mpflag[which(markers=="CD163")]] <- "Myeloid cells"
  d<- cbind(d, coarse_phen_vec)
  ldata[[i]] <- d
}
tsn.marker <- c(immune.markers, "CD31", "SMA", "CD11b", "CD44", "CD31")

#save(ldata, file = "listdatacysift.RData")

#load("/Users/aliamiryousefi/listdatacysift.RData")
phen_vec <- c("Others",  "T cells",  "CD4 T cells",  "Myeloid cells",  "B cells",  "CD8 T cells",  "Endothelial cells",  "Macrophages",  "Tregs",  "M2 Macrophages", "Neutrophils", "Myfibroblasts", "Dendritic cells", "Antigen-presenting cells")
coarse_phen_vec <- c("Others", "T cells", "B cells", "Myeloid cells", "Endothelial cells", "Stromal cells")
lpheprop<- matrix(0, nrow = length(low), ncol = length(coarse_phen_vec))
rownames(lpheprop)<- low
colnames(lpheprop)<- coarse_phen_vec
ic<- 1
for (i in low){
  d<- ldata[[i]]
  den<- length(d$coarse_phen_vec) 
  jc<-1
  for (j in coarse_phen_vec){
    lpheprop[ic, jc]<- sum(d$coarse_phen_vec==j)/ den
    jc<- jc+1
  }
  ic<- ic+1
}

hpheprop<- matrix(0, nrow = length(high), ncol = length(coarse_phen_vec))
rownames(hpheprop)<- high
colnames(hpheprop)<- coarse_phen_vec
ic<- 1
for (i in high){
  d<- ldata[[i]]
  den<- length(d$coarse_phen_vec) 
  jc<-1
  for (j in coarse_phen_vec){
    hpheprop[ic, jc]<- sum(d$coarse_phen_vec==j)/ den
    jc<- jc+1
  }
  ic<- ic+1
}

lpheprop<- as.data.frame(lpheprop)
rownames(lpheprop)<- low
hpheprop<- as.data.frame(hpheprop)
rownames(hpheprop)<- high

lpheprop<- lpheprop[, -c(1,12)]
hpheprop<- hpheprop[, -c(1,12)]



# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
# Combine the data into a single data frame with a 'Group' column
hpheprop_long <- hpheprop %>%
  pivot_longer(cols = everything(), names_to = "Cell_Type", values_to = "Value") %>%
  mutate(Group = "High")

lpheprop_long <- lpheprop %>%
  pivot_longer(cols = everything(), names_to = "Cell_Type", values_to = "Value") %>%
  mutate(Group = "Low")

combined_data <- bind_rows(hpheprop_long, lpheprop_long)

# Calculate means for each cell type and group
mean_values <- combined_data %>%
  group_by(Cell_Type, Group) %>%
  summarise(mean_value = mean(Value), .groups = 'drop')



library(ggplot2)
library(dplyr)

# Define original order of Cell_Type
original_order <- c("T cells", "B cells", "Myeloid cells", "Endothelial cells", "Stromal cells")  # Replace with your actual Cell_Type values

# Convert Cell_Type to a factor with the specified levels
combined_data$Cell_Type <- factor(combined_data$Cell_Type, levels = original_order)

# Plot using ggplot2
ggplot(combined_data, aes(x = Cell_Type, y = Value, fill = Group)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = c("High" = "#7570b3", "Low" = "#1b9e77")) +
  labs(title = "Average Proportions of Cell Types on the Whole Tissue",
       x = "Cell Type", y = "Proportion", fill = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),                         # Increase y-axis text size
    axis.title.x = element_text(size = 16),                        # Increase x-axis title size
    axis.title.y = element_text(size = 16),                        # Increase y-axis title size
    plot.title = element_text(size = 20)                           # Increase title size
  )



###now doing it on the Tumor and on the stroma separetly 


lpheprop<- matrix(0, nrow = length(low), ncol = length(coarse_phen_vec))
rownames(lpheprop)<- low
colnames(lpheprop)<- coarse_phen_vec
ic<- 1
for (i in low){
  d<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8)
  den<- length(d$coarse_phen_vec) 
  jc<-1
  for (j in coarse_phen_vec){
    lpheprop[ic, jc]<- sum(d$coarse_phen_vec==j)/ den
    jc<- jc+1
  }
  ic<- ic+1
}

hpheprop<- matrix(0, nrow = length(high), ncol = length(coarse_phen_vec))
rownames(hpheprop)<- high
colnames(hpheprop)<- coarse_phen_vec
ic<- 1
ssss<-high
high<- high[-10]
for (i in high){
  d<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8)
  den<- length(d$coarse_phen_vec) 
  jc<-1
  for (j in coarse_phen_vec){
    hpheprop[ic, jc]<- sum(d$coarse_phen_vec==j)/ den
    jc<- jc+1
  }
  ic<- ic+1
}


lpheprop<- as.data.frame(lpheprop)
rownames(lpheprop)<- low
hpheprop<- as.data.frame(hpheprop)
rownames(hpheprop)<- ssss

lpheprop<- lpheprop[, -c(1,12)]
hpheprop<- hpheprop[, -c(1,12)]


# Combine the data into a single data frame with a 'Group' column
hpheprop_long <- hpheprop %>%
  pivot_longer(cols = everything(), names_to = "Cell_Type", values_to = "Value") %>%
  mutate(Group = "High")

lpheprop_long <- lpheprop %>%
  pivot_longer(cols = everything(), names_to = "Cell_Type", values_to = "Value") %>%
  mutate(Group = "Low")

combined_data <- bind_rows(hpheprop_long, lpheprop_long)

# Calculate means for each cell type and group
mean_values <- combined_data %>%
  group_by(Cell_Type, Group) %>%
  summarise(mean_value = mean(Value), .groups = 'drop')


# Define original order of Cell_Type
original_order <- c("T cells", "B cells", "Myeloid cells", "Endothelial cells", "Stromal cells")  # Replace with your actual Cell_Type values

# Convert Cell_Type to a factor with the specified levels
combined_data$Cell_Type <- factor(combined_data$Cell_Type, levels = original_order)

# Plot using ggplot2
ggplot(combined_data, aes(x = Cell_Type, y = Value, fill = Group)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = c("High" = "#7570b3", "Low" = "#1b9e77")) +
  labs(title = "Average Proportions of Cell Types on Tumor",
       x = "Cell Type", y = "Proportion", fill = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),                         # Increase y-axis text size
    axis.title.x = element_text(size = 16),                        # Increase x-axis title size
    axis.title.y = element_text(size = 16),                        # Increase y-axis title size
    plot.title = element_text(size = 20)                           # Increase title size
  )
high<-ssss

###now on the stroma


###now doing it on the Tumor and on the stroma separetly 


lpheprop<- matrix(0, nrow = length(low), ncol = length(coarse_phen_vec))
rownames(lpheprop)<- low
colnames(lpheprop)<- coarse_phen_vec
ic<- 1
for (i in low){
  d<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) != 8)
  den<- length(d$coarse_phen_vec) 
  jc<-1
  for (j in coarse_phen_vec){
    lpheprop[ic, jc]<- sum(d$coarse_phen_vec==j)/ den
    jc<- jc+1
  }
  ic<- ic+1
}

hpheprop<- matrix(0, nrow = length(high), ncol = length(coarse_phen_vec))
rownames(hpheprop)<- high
colnames(hpheprop)<- coarse_phen_vec
ic<- 1
for (i in high){
  d<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) != 8)
  den<- length(d$coarse_phen_vec) 
  jc<-1
  for (j in coarse_phen_vec){
    hpheprop[ic, jc]<- sum(d$coarse_phen_vec==j)/ den
    jc<- jc+1
  }
  ic<- ic+1
}

lpheprop<- as.data.frame(lpheprop)
rownames(lpheprop)<- low
hpheprop<- as.data.frame(hpheprop)
rownames(hpheprop)<- high

lpheprop<- lpheprop[, -c(1,12)]
hpheprop<- hpheprop[, -c(1,12)]


# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
# Combine the data into a single data frame with a 'Group' column
hpheprop_long <- hpheprop %>%
  pivot_longer(cols = everything(), names_to = "Cell_Type", values_to = "Value") %>%
  mutate(Group = "High")

lpheprop_long <- lpheprop %>%
  pivot_longer(cols = everything(), names_to = "Cell_Type", values_to = "Value") %>%
  mutate(Group = "Low")

combined_data <- bind_rows(hpheprop_long, lpheprop_long)

# Calculate means for each cell type and group
mean_values <- combined_data %>%
  group_by(Cell_Type, Group) %>%
  summarise(mean_value = mean(Value), .groups = 'drop')


library(ggplot2)
library(dplyr)

# Define original order of Cell_Type
original_order <- c("T cells", "B cells", "Myeloid cells", "Endothelial cells", "Stromal cells")  # Replace with your actual Cell_Type values

# Convert Cell_Type to a factor with the specified levels
combined_data$Cell_Type <- factor(combined_data$Cell_Type, levels = original_order)

# Plot using ggplot2
ggplot(combined_data, aes(x = Cell_Type, y = Value, fill = Group)) +
  stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.9), color = "black") +
  scale_fill_manual(values = c("High" = "#7570b3", "Low" = "#1b9e77")) +
  labs(title = "Average Proportions of Cell Types on Stroma",
       x = "Cell Type", y = "Proportion", fill = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),                         # Increase y-axis text size
    axis.title.x = element_text(size = 16),                        # Increase x-axis title size
    axis.title.y = element_text(size = 16),                        # Increase y-axis title size
    plot.title = element_text(size = 20)                           # Increase title size
  )





#   
#   All Cells
# │
# ├── Immune Cells
# │   ├── CD3d+ (T cells)
# │   │   ├── CD4+ (CD4 T cells)
# │   │   │   └── FOXP3+ (Tregs)
# │   │   └── CD8a+ (CD8 T cells)
# │   │       └── GranzymeB+ (Cytotoxic T cells)
# │   ├── CD19+ / CD20+ (B cells)
# │   ├── CD11b+ / CD68+ / CD11c+ / CD163+ (Myeloid Cells)
# │   │   ├── CD11b+ / CD68+ (Macrophages/DCs)
# │   │   │   ├── CD163+ (M2 Macrophages)
# │   │   │   │   └── CD206+ (Alternative Macrophages)
# │   │   │   └── CD11c+ (Dendritic Cells)
# │   ├── CD15+ (Neutrophils)
# │   └── HLA-DR+ (Antigen-presenting cells)
# │
# ├── CD31+ (Endothelial Cells)
# │
# └── aSMA+ (Myofibroblasts)


#   
#   
#   All Cells
# │
# ├── Epithelial Cells (Prostate Cancer Cells)
# │   ├── HMWCK+ (Basal Cells)
# │   ├── AMCAR+ (Prostate Cancer Cells)
# │   └── Ki67+ (Proliferating Cells)
# │
# ├── Immune Cells
# │   ├── CD3d+ (T cells)
# │   │   ├── CD4+ (CD4 T cells)
# │   │   │   └── FOXP3+ (Tregs)
# │   │   └── CD8a+ (CD8 T cells)
# │   │       └── GranzymeB+ (Cytotoxic T cells)
# │   │       └── PD1+ (Exhausted T cells)
# │   ├── CD19+ / CD20+ (B cells)
# │   ├── CD11b+ / CD68+ / CD11c+ / CD163+ (Myeloid Cells)
# │   │   ├── CD11b+ / CD68+ (Macrophages/DCs)
# │   │   │   ├── CD163+ (M2 Macrophages)
# │   │   │   │   └── CD206+ (Alternative Macrophages)
# │   │   │   └── CD11c+ (Dendritic Cells)
# │   ├── CD15+ (Neutrophils)
# │   └── HLA-DR+ (Antigen-presenting cells)
# │
# ├── Endothelial Cells
# │   └── CD31+ (Endothelial Cells)
# │
# ├── Stromal Cells
# │   ├── aSMA+ (Myofibroblasts)
# │   └── CD44+ (Cancer-associated fibroblasts)
# │
# └── Other Markers
# ├── HLAA+ (MHC class I molecules)
# ├── pTBK1+ (Activated signaling pathways)
# ├── CD103+ (Tissue-resident memory T cells)
# └── TCF1+ (T cell factor 1, associated with T cells)



###

hborder<- matrix(0, nrow = length(high), ncol = length(49:138))
rownames(hborder)<- high
colnames(hborder)<- names(ldata[[1]])[49:138]

ic<- 1
for (i in high){
  l<- as.data.frame(subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 12) == 12))
  jc<-1
  for (j in colnames(hborder)){
    if(strsplit(colnames(hborder)[jc], "_")[[1]][2] == "35r"){
      hborder[ic, jc] <- sum(l[, which(colnames(l)==j)]) / (dim(l)[1] * pi * 35^2)
    }
    if(strsplit(colnames(hborder)[jc], "_")[[1]][2] == "50r"){
      hborder[ic, jc] <- sum(l[, which(colnames(l)==j)]) / (dim(l)[1] * pi * 50^2)
    }
    else if(strsplit(colnames(hborder)[jc], "_")[[1]][2] == "100r"){
      hborder[ic, jc] <- sum(l[, which(colnames(l)==j)]) / (dim(l)[1] * pi * 100^2)
    }
    jc<- jc+1
  }
  ic<- ic+1
  print(hborder)
}



lborder<- matrix(0, nrow = length(low), ncol = length(49:138))
rownames(lborder)<- low
colnames(lborder)<- names(ldata[[1]])[49:138]

ic<- 1
for (i in low){
  l<- as.data.frame(subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 12) == 12))
  jc<-1
  for (j in colnames(lborder)){
    if(strsplit(colnames(lborder)[jc], "_")[[1]][2] == "35r"){
      lborder[ic, jc] <- sum(l[, which(colnames(l)==j)])/ (dim(l)[1] * pi * 35^2)
    }
    if(strsplit(colnames(lborder)[jc], "_")[[1]][2] == "50r"){
      lborder[ic, jc] <- sum(l[, which(colnames(l)==j)]) / (dim(l)[1] * pi * 50^2)
    }
    if(strsplit(colnames(lborder)[jc], "_")[[1]][2] == "100r"){
      lborder[ic, jc] <- sum(l[, which(colnames(l)==j)]) / (dim(l)[1] * pi * 100^2)
    }
    jc<- jc+1
  }
  ic<- ic+1
  print(lborder)
}

lborder<- lborder[-c(8,11,12),]

lborder<- t(lborder)
hborder<- t(hborder)



par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
for (i in rownames(hborder)[1:21]){
  # Calculate the percentages
  l_percent <- 100 * as.numeric(lborder[which(rownames(lborder)== i), ])
  h_percent <- 100 * as.numeric(hborder[which(rownames(lborder)== i), ])
  
  # Combine the data
  combined_data <- c(l_percent, h_percent)
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 10), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#7570b3", "#1b9e77"), outline=FALSE)
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups ,              # Data
             method = "jitter", # Random noise
             pch = 19,          # Pch symbols
             col = 1,           # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE)        # Add it over
  
  t_test_result <- wilcox.test(l_percent, h_percent, alternative = "two.sided", exact = FALSE, correct = FALSE)
  
  # Add t-test result as subtitle
  p_value <- t_test_result$p.value
  if(p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }
}


dev.off()

par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
for (i in rownames(hborder)[22:42]){
  # Calculate the percentages
  l_percent <- 100 * as.numeric(lborder[which(rownames(lborder)== i), ])
  h_percent <- 100 * as.numeric(hborder[which(rownames(lborder)== i), ])
  
  # Combine the data
  combined_data <- c(l_percent, h_percent)
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 10), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#7570b3", "#1b9e77"), outline=FALSE)
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups ,              # Data
             method = "jitter", # Random noise
             pch = 19,          # Pch symbols
             col = 1,           # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE)        # Add it over
  
  t_test_result <- wilcox.test(l_percent, h_percent, alternative = "two.sided", exact = FALSE, correct = FALSE)
  
  # Add t-test result as subtitle
  p_value <- t_test_result$p.value
  if(p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }
}

par(mfrow = c(3, 7), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
for (i in rownames(hborder)[43:63]){
  # Calculate the percentages
  l_percent <- 100 * as.numeric(lborder[which(rownames(lborder)== i), ])
  h_percent <- 100 * as.numeric(hborder[which(rownames(lborder)== i), ])
  
  # Combine the data
  combined_data <- c(l_percent, h_percent)
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 10), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#7570b3", "#1b9e77"), outline=FALSE)
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups ,              # Data
             method = "jitter", # Random noise
             pch = 19,          # Pch symbols
             col = 1,           # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE)        # Add it over
  
  t_test_result <- wilcox.test(l_percent, h_percent, alternative = "two.sided", exact = FALSE, correct = FALSE)
  
  # Add t-test result as subtitle
  p_value <- t_test_result$p.value
  if(p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }
}
par(mfrow = c(4, 7), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
for (i in rownames(hborder)[64:90]){
  # Calculate the percentages
  l_percent <- 100 * as.numeric(lborder[which(rownames(lborder)== i), ])
  h_percent <- 100 * as.numeric(hborder[which(rownames(lborder)== i), ])
  
  # Combine the data
  combined_data <- c(l_percent, h_percent)
  
  # Create a vector to represent the groups (L and H)
  groups <- c(rep("L", 10), rep("H", 13))
  
  # Create boxplot with colored boxes and percentage symbol on the y-axis
  boxplot(combined_data ~ groups, range=3,
          main = paste(i, "+", sep = ""), 
          ylab = "", xlab = "", yaxt = "n", 
          col = c("#7570b3", "#1b9e77"), outline=FALSE)
  axis(2, at = pretty(range(combined_data)), labels = paste0(pretty(range(combined_data)), "%"), las=2)
  stripchart(combined_data ~ groups ,              # Data
             method = "jitter", # Random noise
             pch = 19,          # Pch symbols
             col = 1,           # Color of the symbol
             vertical = TRUE,   # Vertical mode
             add = TRUE)        # Add it over
  
  t_test_result <- wilcox.test(l_percent, h_percent, alternative = "two.sided", exact = FALSE, correct = FALSE)
  
  # Add t-test result as subtitle
  p_value <- t_test_result$p.value
  if(p_value < 0.05) {
    subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
  } else {
    subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
  }
}



##### DBSCAN of the immune cells on the tumor
# Load libraries

#trying it for one of the images

T.inf<- numeric(length(names(ldata)))
for (i in names(ldata)){
  tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8 )
  if (dim(tumor)[1] > 50){
    data<- subset(tumor, tumor$phen_vec=="T cells")
    ppp <- as.ppp(data[, c("x", "y")], W = owin(c(min(data$x), max(data$x)), c(min(data$y), max(data$y))))
    K <- Lest(ppp)
    T.inf[which(names(ldata)==i)] <-  mad.test(ppp, Lest, nsim= 5, rmax=10, use.theo=TRUE, interpolate = TRUE)$p.value
    print(T.inf[which(names(ldata)==i)])
    plot(K, main = "Ripley's K Function")
    #str(ppp)
   # Sys.sleep(2)
  }
}

###go with the 200/500*200/500 wondow size and 10 samples of the tumor are
pdf("k-function_1k_windowsize-envelope.pdf", width=12.5, height=8.5)
ws<-2000
T.inf<- numeric(length(names(ldata)))
par(mfrow = c(4, 4), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
for (i in names(ldata)[1:2]){
  tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8 & ldata[[i]]$tls_id == 0 )
  of.tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) != 8 & ldata[[i]]$tls_id == 0)
  if (dim(tumor)[1] > 50){
    c<-1
    data<- subset(tumor, tumor$phen_vec=="T cells")
    for (j in 1: 20){
      xstart<- sample(c(min(tumor$x):(max(tumor$x)-ws)), 1)
      ystart<-  sample(c(min(tumor$y):(max(tumor$y)-ws)), 1)
      datar<- subset(tumor, tumor$x < xstart+ ws & tumor$x > xstart & tumor$y < ystart + ws & tumor$y > ystart)
      if (dim(datar)[1] > 50 && c < 5){
        c<- c+1
        ppp <- as.ppp(datar[, c("x", "y")], W = owin(c(min(datar$x), max(datar$x)), c(min(datar$y), max(datar$y))))
        #K <- Lest(ppp, rmax=300/0.325)
        #T.inf[which(names(ldata)==i)] <-  mad.test(ppp, Lest, nsim= 10, rmax=10, use.theo=TRUE, interpolate = TRUE)$p.value
        #print(T.inf[which(names(ldata)==i)])
        plot(envelope(ppp, Lest, nsim= 5, fix.n=TRUE, global = TRUE, rmax= 150/0.325), main = paste(paste(paste(paste(i, "Tumor at x:", sep=" "), xstart, sep=""), "y:", sep=""), ystart, sep=""))
      }
      else {
      }
    }
  }
      if (dim(of.tumor)[1] > 50){
        c<-1
        data<- subset(of.tumor, of.tumor$phen_vec=="T cells")
        #on the stroma
        for (j in 1: 20){
        xstart<- sample(c(min(of.tumor$x):(max(of.tumor$x)-ws)), 1)
        ystart<-  sample(c(min(of.tumor$y):(max(of.tumor$y)-ws)), 1)
        datar<- subset(of.tumor, of.tumor$x < xstart+ ws & of.tumor$x > xstart & of.tumor$y < ystart + ws & of.tumor$y > ystart)
        if (dim(datar)[1] > 50 && c < 5){
          c<- c+1
          ppp <- as.ppp(datar[, c("x", "y")], W = owin(c(min(datar$x), max(datar$x)), c(min(datar$y), max(datar$y))))
          #K <- Lest(ppp, rmax= 300/0.325)
          #T.inf[which(names(ldata)==i)] <-  mad.test(ppp, Lest, nsim= 10, rmax=10, use.theo=TRUE, interpolate = TRUE)$p.value
          # print(T.inf[which(names(ldata)==i)])
          plot(envelope(ppp, Lest, nsim= 5, fix.n=TRUE, global = TRUE, rmax=150/0.325), main = paste(paste(paste(paste(i, "Stroma at x:", sep=" "), xstart, sep=""), "y;", sep=""), ystart, sep=""))
    # Sys.sleep(2)
        }
        else {}
    }  
  }
}
dev.off()
dev.off()


#####trying it on the non-window and instead on the whole tumor area
mean.L<- numeric(length(names(ldata)))
pdf("k-function_1k_windowsize-envelope.pdf", width=12.5, height=8.5)
ws<-150/0.325
T.inf<- numeric(length(names(ldata)))
par(mfrow = c(4, 4), mar = c(2.2, 4.1, 4.1, 2.1)) #mar = c(bottom, left, top, right)  the default is c(5.1, 4.1, 4.1, 2.1)
mean.holder <- matrix(0, 10, length(names(ldata)))
colnames(mean.holder)<- names(ldata)
for (i in names(ldata)){
  tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) == 8 & ldata[[i]]$tls_id == 0 )
  of.tumor<- subset(ldata[[i]], bitAnd(ldata[[i]]$cflag, 8) != 8 & ldata[[i]]$tls_id == 0)
  mean.holderv<- numeric(10)
  if (dim(tumor)[1] > 50){
    c<-1
    data<- subset(tumor, tumor$phen_vec=="T cells")
    for (j in 1: 50){
      xstart<- sample(c(min(tumor$x):(max(tumor$x)-ws)), 1)
      ystart<-  sample(c(min(tumor$y):(max(tumor$y)-ws)), 1)
      datar<- subset(tumor, tumor$x < xstart+ ws & tumor$x > xstart & tumor$y < ystart + ws & tumor$y > ystart)
      if (dim(datar)[1] > 50 && c < 11){
        c<- c+1
        ppp <- as.ppp(datar[, c("x", "y")], W = owin(c(min(data$x), max(data$x)), c(min(data$y), max(data$y))))
        L <- Lest(ppp, rmax=150/0.325)
        mean.holder[(c-1), which(names(ldata)==i)] <- mean(L$border-L$theo)
        mean.holderv[(c-1)] <- mean(L$border-L$theo)
        #print(head(L$theo))
        #K <- Lest(ppp, rmax=300/0.325)
        #T.inf[which(names(ldata)==i)] <-  mad.test(ppp, Lest, nsim= 10, rmax=10, use.theo=TRUE, interpolate = TRUE)$p.value
        #print(T.inf[which(names(ldata)==i)])
        #plot(envelope(ppp, Lest, nsim= 5, fix.n=TRUE, global = TRUE, rmax= 150/0.325), main = paste(paste(paste(paste(i, "Tumor at x:", sep=" "), xstart, sep=""), "y:", sep=""), ystart, sep=""))
      }
      else {
      }
    }
  }
  
  mean.L[which(names(ldata)==i)]<- mean(mean.holderv[mean.holderv != 0 & !is.na(mean.holderv)], na.rm = TRUE)
  print(mean.holder)
  #if (dim(of.tumor)[1] > 50){
  #  c<-1
  #  data<- subset(of.tumor, of.tumor$phen_vec=="T cells")
  #  #on the stroma
  #  for (j in 1: 20){
  #   xstart<- sample(c(min(of.tumor$x):(max(of.tumor$x)-ws)), 1)
  #   ystart<-  sample(c(min(of.tumor$y):(max(of.tumor$y)-ws)), 1)
  #   datar<- subset(of.tumor, of.tumor$x < xstart+ ws & of.tumor$x > xstart & of.tumor$y < ystart + ws & of.tumor$y > ystart)
  #    if (dim(datar)[1] > 50 && c < 5){
  #      c<- c+1
  #      ppp <- as.ppp(datar[, c("x", "y")], W = owin(c(min(data$x), max(data$x)), c(min(data$y), max(data$y))))
  #     #K <- Lest(ppp, rmax= 300/0.325)
  #      #T.inf[which(names(ldata)==i)] <-  mad.test(ppp, Lest, nsim= 10, rmax=10, use.theo=TRUE, interpolate = TRUE)$p.value
  #      # print(T.inf[which(names(ldata)==i)])
  #      plot(envelope(ppp, Lest, nsim= 5, fix.n=TRUE, global = TRUE, rmax=150/0.325), main = paste(paste(paste(paste(i, "Stroma at x:", sep=" "), xstart, sep=""), "y;", sep=""), ystart, sep=""))
  #      # Sys.sleep(2)
  #    }
  #    else {}
  #  }  
  #}
}

dev.off()
dev.off()


names(mean.L)<- names(ldata)

colors <- ifelse(names(sort(mean.L, decreasing = TRUE)) %in% high, "#7570b3", "#1b9e77")
barplot(sort(mean.L, decreasing = TRUE), las = 2, col = colors, main = "Average 10 Reply's L clustering on 200x200 mu", ylab = "", xlab = "")
# Adding a legend
legend("topright", legend = c("High", "Low"), fill = c("#7570b3", "#1b9e77"), title = "Categories")

data <- data.frame(
  value = mean.L,
  group = ifelse(names(mean.L) %in% high, "High", "Low"),
  name = names(mean.L)
)

# Perform t-test
t_test_result <- t.test(value ~ group, data = data)

# Extract p-value
p_value <- t_test_result$p.value

# Create the box plot with overlaid points
p <- ggplot(data, aes(x = group, y = value)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA) +    # Box plot without outliers
  geom_jitter(width = 0.2, aes(color = group), size = 3, shape = 21, fill = NA, stroke = 1, color = "black") + # Outline points
  #geom_jitter(width = 0.2, aes(color = group), size = 3) + # Color points
  scale_fill_manual(values = c("High" = "#7570b3", "Low" = "#1b9e77")) + # Custom fill colors
  scale_color_manual(values = c("High" = "#7570b3", "Low" = "#1b9e77")) + # Custom point colors
  labs(title = "Reply's L test",
       x = "Group",
       y = "Value") +
  theme_minimal() +
  geom_text(aes(x = 1.5, y = 3600, 
                label = paste("p-value =", format(p_value, digits = 3))),
            color = ifelse(p_value < 0.05, "red", "black"), size = 5)

# Print the plot
print(p)


#### project pursuit

library(fastICA)
library(umap)
if (!requireNamespace("uwot", quietly = TRUE)) {
  install.packages("uwot")
}
library(uwot)
library(ggplot2)
###tSNE or leiden clustering with plotting dimensions
tsn.marker <- c(immune.markers, "CD31", "SMA", "CD11b", "CD44", "CD31", "coarse_phen_vec")

 for (i in 1:length(names(ldata))){
   data<- as.data.frame(ldata[[i]])
   inp<- data[c(which(data$coarse_phen_vec!="Others")), c(which(colnames(data)%in%tsn.marker))] 
   inp1<- inp[c(which(inp$coarse_phen_vec!="Stromal cells")),]
   inp2<- inp[c(which(inp$coarse_phen_vec=="Stromal cells")),]
   inp2<- inp2[sample(dim(inp2)[1], round(dim(inp2)[1]/3)), ]
   #inp2<- inp2[1:round(dim(inp2)[1]/3), ]
   inp<- rbind(inp1, inp2)
   
   #ica<- fastICA(inp[, -which(colnames(inp)=="coarse_phen_vec")], 5)
   
   #plot(ica$S, mainb="ICA components")
   
   
   system.time({
     umap_result <- umap(inp[, -which(colnames(inp)=="coarse_phen_vec")], n_neighbors = 8, min_dist = 0.5, n_threads = 8)
   })
   
   umap_df <- data.frame(UMAP1 = umap_result[, 1], UMAP2 = umap_result[, 2], Color = inp$coarse_phen_vec)
   
   # Plot with ggplot2
   p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = factor(Color))) +
     geom_point(size = 0.2, alpha = 0.5) +  # Adjust size and alpha for smaller, more transparent points
     theme_minimal() +
     labs(title = "UMAP Clustering", x = "UMAP1", y = "UMAP2", color = "Cluster") +
     guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
   
   p
   
   file_name <- paste0("UMAP_Plot_", names(ldata)[i], ".pdf")
   
   # Save the plot
   ggsave(file_name, plot = p, width = 8, height = 6)
 }





LSP33<- ldata[["LSP12633"]]
LSP35<- ldata[["LSP12635"]]
tumor<- subset(LSP33, bitAnd(LSP33$cflag, 8) == 8 )

immune.markers

library(spatstat)

# Load your data (assuming it's in a CSV file)
data <- tumor

# Step 3: Spatial Statistics

data<- subset(data, data$coarse_phen_vec=="B cells")
ppp <- as.ppp(data[, c("x", "y")], W = owin(c(min(data$x), max(data$x)), c(min(data$y), max(data$y))))
K <- Kest(ppp)
E<- envelope(ppp, Kest, nsim= 100, fix.n=TRUE, global = TRUE)
#nongraphical tests
#mad.test(ppp, Kest, nsim= 100, rmax=100, use.theo=TRUE)
dclf.test(ppp, Kest, nsim= 100, rmax=10, use.theo=TRUE)$p.value
plot(K, main = "Ripley's K Function")
plot(E)


##### the KiGate analysis

#####With the latest ldata with the Ki67 column

sums <- lapply(ldata, function(df) sum(df[["KiGate"]])/sum(df[["tls_id"]]!=0))
barplot(sort(unlist(sums), decreasing = TRUE), las=2)
sorted_sums <- sort(unlist(sums), decreasing = TRUE)
colors <- ifelse(names(sorted_sums) %in% high, "#7570b3", "#1b9e77")
barplot(sorted_sums, las = 2, col = colors, main = "Proportion of Ki67+ in TLSs", ylab = "", xlab = "")
# Adding a legend
legend("topright", legend = c("High", "Low"), fill = c("#7570b3", "#1b9e77"), title = "Categories")


low_values <- na.omit(sums[names(sums) %in% low])
low_values<-subset(low_values, low_values>0)
high_values <- na.omit(sums[names(sums) %in% high])
t_test_result <- t.test(unlist(low_values), unlist(high_values))
print(t_test_result)
data <- data.frame(value = c(unlist(low_values), unlist(high_values)),group = factor(rep(c("Low", "High"), c(length(low_values), length(high_values)))))

library(ggplot2)
ggplot(data, aes(x = group, y = value, color = group)) +
  geom_boxplot(outlier.shape = NA) +  # Avoid double plotting of outliers
  geom_jitter(width = 0.2, size = 2) +  # Add jittered points
  theme_minimal() +
  labs(title = "Boxplots with Overlaid Data Points",
       x = "Group", y = "Values") +
  scale_color_manual(values = c("Low" = "#1b9e77", "High" = "#7570b3")) +
  theme(legend.position = "none")



h.ntls<- numeric(length(high))
c<- 1
for (i in high){l<- ldata[[which(names(ldata)== i )]]
h.ntls[c]<- length(unique(l$tls_id))
c<- c+1
}
h.ntls<- h.ntls-1

l.ntls<- numeric(length(low))
c<- 1
for (i in low){l<- ldata[[which(names(ldata)== i )]]
l.ntls[c]<- length(unique(l$tls_id))
c<- c+1
}
l.ntls<- l.ntls-1


markers<- "Ki67"

h.mat<- as.data.frame(matrix(0, length(markers), sum(h.ntls)))
#colnames(h.mat)<- high
rownames(h.mat)<- markers
c<-1
aa<-1
for (j in high){l<- ldata[[which(names(ldata)== j )]]
for (t in 1:h.ntls[c]){
  l<-subset(l, l$tls_id == c(unique(l$tls_id)[-which(unique(l$tls_id)==0)])[t])
  if (dim(l)[1] > 0){
    h.markprop<- numeric(length(markers))
    for (i in markers) {
      sum<-0
      for (k in 1:dim(l)[1]){
        sum<- l$KiGate[k] + sum
      }
      h.markprop[which(markers == i)] <- sum/dim(l)[1]
    }
    h.mat[,aa]<- h.markprop
  }
  aa<- aa + 1
  print(h.mat)
  l<- ldata[[which(names(ldata)== j )]]
}
c<- c+1
}

l.mat<- as.data.frame(matrix(0, length(markers), sum(l.ntls)))
#colnames(h.mat)<- high
rownames(l.mat)<- markers
c<-1
aa<-1
for (j in low){l<- ldata[[which(names(ldata)== j )]]
for (t in 1:l.ntls[c]){
  l<-subset(l, l$tls_id == c(unique(l$tls_id)[-which(unique(l$tls_id)==0)])[t])
  if (dim(l)[1] > 0){
    l.markprop<- numeric(length(markers))
    for (i in markers) {
      sum<-0
      for (k in 1:dim(l)[1]){
        sum<- l$KiGate[k] + sum
      }
      l.markprop[which(markers == i)] <- sum/dim(l)[1]
    }
    l.mat[,aa]<- l.markprop
  }
  aa<- aa + 1
  print(l.mat)
  l<- ldata[[which(names(ldata)== j )]]
}
c<- c+1
}
ll.mat<- l.mat
l.mat<- as.data.frame(t(subset(t(l.mat), t(l.mat) < 0.4)))

pdf("Ki67.pdf")
# Calculate the percentages
l_percent <- 100 * as.numeric(l.mat[which(rownames(l.mat)== i), ])
h_percent <- 100 * as.numeric(h.mat[which(rownames(l.mat)== i), ])

# Combine the data
combined_data <- c(l_percent, h_percent)

# Create a vector to represent the groups (L and H)
groups <- c(rep("L", length(l.mat)), rep("H", 190))

# Create boxplot with colored boxes and percentage symbol on the y-axis
boxplot(combined_data ~ groups, range=4,
        main = paste(i, "+", sep = ""), 
        ylab = "", xlab = "", yaxt = "n", 
        col = c("#1b9e77", "#7570b3"), outline=FALSE)

# Adjust the y-axis labels to ensure at least two ticks are printed
axis(2, at = pretty(range(combined_data), n = round(exp(1/median(combined_data)))+10), labels = paste0(pretty(range(combined_data), n = round(exp(1/median(combined_data)))+10), "%"), las=2)
#print(mean(combined_data))
# Generate distinctive colors for each group
colors <- c(rainbow(26)[1:13], rainbow(26)[14:26])

# Replicate colors for each group
group_colors <- c(rep(colors[1:13], l.ntls), rep(colors[14:26], h.ntls))

# Plot each data point individually with jitteriness and respective color
for (j in 1:length(combined_data)) {
  jittered_x <- ifelse(groups[j] == "H", 1, 2) + runif(1, -0.2, 0.2)  # Add jitter to x-coordinate
  points(jittered_x, combined_data[j], pch = 20, col = adjustcolor(group_colors[j], alpha.f = 0.5))
}

# Perform t-test
t_test_result <- t.test(l_percent, h_percent)

# Add t-test result as subtitle
p_value <- t_test_result$p.value
if (p_value < 0.05) {
  subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
  mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
} else {
  subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
  mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
}
l.mat<- ll.mat
dev.off()




##########exaustive search for aggregates

d<- ldata[["LSP12657"]]


library(spatstat) # Ensure you have the necessary library for as.ppp, owin, Lest functions

imm.cruse <- function(ws, LSP, index, phen, fig = TRUE, creep = FALSE) {
  # Window size (ws) in micron, LSP identifies in the list of ldata object,
  # index to be calculated, phenotype for the indexing e.g B cells (reading off the coarse_phen_vec)
  L.models<- list()
  pws <- ws / 0.325 # pixel window size
  d <- ldata[[LSP]]
  xstart <- 1
  ystart <- 1
  
  if (fig==TRUE){
    a<- max(d$x); b<- max(d$y)
    par(mar = c((10.66-(b/1500)), (10.66-(a/1500)), (10.66-(b/1500)), (10.66-(a/1500))))
    plot(d$x, d$y, pch=19, cex=0.01, cex.axis=1,  # Axis tick labels larger
         cex.lab=1,   # Axis title labels larger
         cex.main=1.4,    # Main title larger 
         col="lightgrey", main = LSP, xlab = paste("Window Size", ws, sep = "="), ylab = "", ylim = range(d$y), xlim =range(d$x), col.main = "navy"
         #,panel.first = grid()
    )  
    if (sum(bitAnd(d$cflag, 8) == 8) > 0){
      ss<- subset(d, bitAnd(d$cflag, 8) == 8)
      #s1<-sample(dim(ss)[1], round(dim(ss)[1]/50))
      #ss<- ss[s1, ]
      points(ss$x, ss$y, col=adjustcolor("lightblue", alpha.f = 0.2), pch=19, cex=0.005)
    }
    if (sum(d$coarse_phen_vec == phen) > 0){
      if (phen == "T cells"){co<- "green"}
      if (phen == "B cells"){co<- "red"}
      ss<- subset(d, d$coarse_phen_vec == phen)
      #s1<-sample(dim(ss)[1], round(dim(ss)[1]/4))
      #ss<- ss[s1, ]
      points(ss$x, ss$y, col=adjustcolor(co, alpha.f = 0.4), pch=19, cex=0.005)
    }
    for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
    for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i + pws/2, col = adjustcolor("darkolivegreen", alpha.f = 0.5), lty = 3, lwd = 1)}
    for (j in 0:ceiling(max(d$x) / pws)){abline(v = xstart + pws * j, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
  }
  
  ny <- ceiling(max(d$y) / pws)
  nx <- ceiling(max(d$x) / pws)
  index.holder <- matrix(0, ny, nx)
  xstart <- 1
  ystart <- 1
  
  # Initialize progress bar
  total_steps <- nx * ny * creep^2
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  step <- 0
  c<-1
  for(k in 1:creep){
    xstart<- k*round(pws/creep)-round(pws/creep)
    for (l in 1:creep){
    ystart<- l*round(pws/creep)-round(pws/creep)
  for (i in 1:nx) {
    for (j in 1:ny) {
      data <- subset(d, d$x < xstart + pws * i & d$x > xstart + pws * (i - 1) & d$y < ystart + pws * j & d$y > ystart + pws * (j - 1))
      if (nrow(data) > 50) { # if the box has more than 10 cells
        phen.data <- subset(data, data$coarse_phen_vec == phen)
        if (nrow(phen.data) > 5) {
          if (index == "K-integral") {
            ppp <- as.ppp(phen.data[, c("x", "y")], W = owin(c(min(phen.data$x), max(phen.data$x)), c(min(phen.data$y), max(phen.data$y))))
            L <- Lest(ppp, rmax = ws)
            index.holder[j, i] <- mean(L$border - L$theo, na.rm = TRUE)
            differences <- L$border - L$theo
            L.models[[c]]<- L
            if(fig){
              if(k + l < 3){
                smoothed <- loess(differences ~ seq_along(differences), span = 0.3)
                smoothed_values <- predict(smoothed, seq_along(differences))
                lines(seq(xstart + pws * (i-1), xstart + pws * (i), pws/513)[-1], c(smoothed_values + pws/2 + pws*(j-1)), col='plum1', lwd = 2*ws/500)
              }
              D<- L$border - L$theo
              Dp<- D[D>0]
              text(x = xstart + pws * (i-0.5), y = ystart + pws/2 + pws*(j-1), labels = round(mean(Dp, na.rm = TRUE)), col = "plum4", cex = ws/700 + ws/700 * (round(mean(Dp, na.rm = TRUE))/700), font = 3)
            # Add the smoothed line to the plot
            #lines(smoothed_values, col = "red", lwd = 2)
            }
            c<- c+1
          }
          
          if (index == "Density") {
            # Add Density calculation here
          }
          if (index == "Entropy") {
            # Add Entropy calculation here
          }
        }
      }
      # Update progress bar
      step <- step + 1
      setTxtProgressBar(pb, step)
    }
    }
    }
  }
  
  # Close progress bar
  close(pb)
  
  return(L.models)
}

immB.list<- list()
immT.list<- list()
for (i in 1:length(names(ldata))){
  print(names(ldata)[i])
  immB.list[[i]] <- imm.cruse(500, names(ldata)[i], "K-integral", "B cells", fig = FALSE, creep = 2)
  immT.list[[i]] <- imm.cruse(500, names(ldata)[i], "K-integral", "T cells", fig = FALSE, creep = 2)
}

for (i in 1:length(immB.list)){
  for (j in 1:length(immB.list[[i]])){
    my_vector <- immB.list[[i]][[j]]$border-immB.list[[i]][[j]]$theo
    
    chunk_size <- 10
    n_chunks <- ceiling(length(my_vector) / chunk_size)
    
    immB.list[[i]][[j]] <- sapply(1:n_chunks, function(i) {
      # Define the start and end indices for the current chunk
      start_index <- (i - 1) * chunk_size + 1
      end_index <- min(i * chunk_size, length(my_vector))
      
      # Extract the current chunk
      current_chunk <- my_vector[start_index:end_index]
      
      # Calculate the mean, ignoring NA values
      mean(current_chunk, na.rm = TRUE)
    })
  }
}

for (i in 1:length(immT.list)){
  for (j in 1:length(immT.list[[i]])){
    my_vector <- immT.list[[i]][[j]]$border-immT.list[[i]][[j]]$theo
    
    chunk_size <- 10
    n_chunks <- ceiling(length(my_vector) / chunk_size)
    
    immT.list[[i]][[j]] <- sapply(1:n_chunks, function(i) {
      # Define the start and end indices for the current chunk
      start_index <- (i - 1) * chunk_size + 1
      end_index <- min(i * chunk_size, length(my_vector))
      
      # Extract the current chunk
      current_chunk <- my_vector[start_index:end_index]
      
      # Calculate the mean, ignoring NA values
      mean(current_chunk, na.rm = TRUE)
    })
  }
}



# Access the specific sublist


immB50<- for (i in 1:length(names(ldata))){
  sublist <- immB.list[[i]]
  
  # Ensure all vectors are of the same length
  vector_length <- min(sapply(sublist, length))
  
  # Truncate vectors to the smallest length (if needed)
  truncated_sublist <- lapply(sublist, function(vec) vec[1:vector_length])
  
  # Calculate the element-wise mean across all vectors
  element_wise_means <- rowMeans(do.call(cbind, truncated_sublist), na.rm = TRUE)
  
  # Print the element-wise means
  print(element_wise_means)
  ldata[[i]]$B.CR500 <- numeric(dim(ldata[[i]])[1])
  ldata[[i]]$B.CR500[1:52]<-  element_wise_means
}


immT50<- for (i in 1:length(names(ldata))){
  sublist <- immT.list[[i]]
  
  # Ensure all vectors are of the same length
  vector_length <- min(sapply(sublist, length))
  
  # Truncate vectors to the smallest length (if needed)
  truncated_sublist <- lapply(sublist, function(vec) vec[1:vector_length])
  
  # Calculate the element-wise mean across all vectors
  element_wise_means <- rowMeans(do.call(cbind, truncated_sublist), na.rm = TRUE)
  
  # Print the element-wise means
  print(element_wise_means)
  ldata[[i]]$T.CR500 <- numeric(dim(ldata[[i]])[1])
  ldata[[i]]$T.CR500[1:52]<-  element_wise_means
}


save(ldata, file = "listdatacysift2_Phen_Ki_BT.RData")
#and load it after that.

##saved the 500 micron Clustering and Regularity values with the 10 micron interval in the first 52 elements of the B and T .CR500 vectors on each dataset respectively!

load("/Users/aliamiryousefi/listdatacysift2_Phen_Ki_BT.RData")




for (i in 1:length(names(ldata))){
  
pdf(paste(names(ldata)[i], ".pdf", sep ="CR500"), width=10.5, height=10.5)
par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(1,1,1,1))

imm.cruse(700, LSP = names(ldata)[i], index = "K-integral", phen = "B cells", fig = TRUE, creep = 1)

par(mar=c(5.1, 4.1, 4.1, 2.1))
example_data<- ldata[[i]]$B.CR500[1:50]
example_data1<-example_data[example_data > 0]

total_sum <- sum(example_data1)
cumulative_sum <- cumsum(example_data1)
half_sum <- total_sum / 2
closest_index <- which.min(abs(cumulative_sum - half_sum))
colors <- ifelse(example_data < 0, "darkgrey", "red")
colors[1:closest_index] <- "darkred"


# Define custom x-axis labels
x_labels <- seq(10, 500, by=10)

# Create the barplot
barplot(
  example_data,
  ylim = c(-100, 250),
  col = colors,
  names.arg = x_labels,
  xlab = "Radius (in micron)",
  ylab = "Clustering index",
  main = paste("B cells C&R500", names(ldata)[i], sep=" for "), las = 2, 
  border = "black" # Add borders to bars
)

# Add horizontal grid lines
abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
abline(h = mean(example_data[which(example_data > 0)]), col = "sienna", lty = 1, lwd = 1.5)
text(x= 45, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("M+", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="sienna3", font= 2)
text(x= closest_index + 3, y = -50,las = 2,  labels = paste("R50", closest_index*10 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="darkred", font= 3, srt = 90)

# Add a legend
legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "red"), border = "black")

imm.cruse(700, LSP = names(ldata)[i], index = "K-integral", phen = "T cells", fig = TRUE, creep = 1)
par(mar=c(5.1, 4.1, 4.1, 2.1))

example_data<- ldata[[i]]$T.CR500[1:50]
example_data1<-example_data[example_data > 0]
total_sum <- sum(example_data1)
cumulative_sum <- cumsum(example_data1)
half_sum <- total_sum / 2
closest_index <- which.min(abs(cumulative_sum - half_sum))
colors <- ifelse(example_data < 0, "darkgrey", "green")
colors[1:closest_index] <- "darkgreen"

# Define custom x-axis labels
x_labels <- seq(10, 500, by=10)

# Create the barplot
barplot(
  example_data,
  ylim = c(-100, 250),
  col = colors,
  names.arg = x_labels,
  xlab = "Radius (in micron)",
  ylab = "Clustering index",
  main = paste("T cells C&R500", names(ldata)[i], sep=" for "), las = 2, 
  border = "black" # Add borders to bars
)

# Add horizontal grid lines
abline(h = seq(-100, 250, by=25), col = "darkgrey", lty = "dotted")
abline(h = mean(example_data[which(example_data > 0)]), col = "olivedrab", lty = 1, lwd = 1.5)
text(x= 45, y = mean(example_data[which(example_data > 0)]) + 10, labels = paste("M+", round(mean(example_data[which(example_data > 0)])), sep = ": "), cex = 1.2, col="olivedrab4", font= 2)
text(x= closest_index + 3, y = -50,las = 2,  labels = paste("R50", closest_index*10 + sample(-5:4, 1), sep = ": "), cex = 1.2, col="darkgreen", font= 3, srt = 90)

#text(, y_pos_outer1, paste(7, "um", " "), cex = 0.7, col="darkgrey", font= 2)
# Add a legend
legend("topright", legend = c("Regular", "Clustered"), fill = c("darkgrey", "green"), border = "black")

dev.off()
}


###saving the M+ and R50 for B and T cells for each of the samples then plot them in a bar plot!
B.low_high.M <- numeric(26)
T.low_high.M <- numeric(26)

B.low_high.R <- numeric(26)
T.low_high.R <- numeric(26)
c<-1
for (i in c(low,high)){
  example_data<- ldata[[i]]$B.CR500[1:50]
  example_data1<-example_data[example_data > 0]
  total_sum <- sum(example_data1)
  cumulative_sum <- cumsum(example_data1)
  half_sum <- total_sum / 2
  closest_index <- which.min(abs(cumulative_sum - half_sum))
  B.low_high.M[c] <- mean(example_data[which(example_data > 0)])
  B.low_high.R[c] <- closest_index*10
  
  example_data<- ldata[[i]]$T.CR500[1:50]
  example_data1<-example_data[example_data > 0]
  total_sum <- sum(example_data1)
  cumulative_sum <- cumsum(example_data1)
  half_sum <- total_sum / 2
  closest_index <- which.min(abs(cumulative_sum - half_sum))
  T.low_high.M[c] <- mean(example_data[which(example_data > 0)])
  T.low_high.R[c] <- closest_index*10

  c<- c+1
}

library(ggplot2)
library(ggpubr)

install.packages("gridExtra")
library(gridExtra)

# Define a function to get the y-axis limit
extend_ylim <- function(values, factor = 1.1) {
  range <- range(values, na.rm = TRUE)
  c(range[1], range[2] * factor)
}

# Create the first plot for B cells, Average clustering index
plot1 <- ggplot(data.frame(Values = B.low_high.M, Condition = factor(rep(c("Low", "High"), c(13, 13)))),
                aes(x = Condition, y = Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 3, color = "red") +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  labs(title = "Clustering in 0.5mm for B cells",
       x = "Gleason", y = "Average deviation") +
  ylim(extend_ylim(B.low_high.M)) +
  theme_minimal(base_size = 15) +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  stat_compare_means(method = "t.test", label = "p.format", 
                     size = 6, vjust = 2)  # Adjust vertical position of the p-value

# Create the second plot for T cells, Average clustering index
plot2 <- ggplot(data.frame(Values = T.low_high.M, Condition = factor(rep(c("Low", "High"), c(13, 13)))),
                aes(x = Condition, y = Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 3, color = "green") +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  labs(title = "Clustering in 0.5mm for T cells",
       x = "Gleason", y = "") +
  ylim(extend_ylim(T.low_high.M)) +
  theme_minimal(base_size = 15) +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  stat_compare_means(method = "t.test", label = "p.format", 
                     size = 6, vjust = 2)  # Adjust vertical position of the p-value

# Create the third plot for B cells, R 50
plot3 <- ggplot(data.frame(Values = B.low_high.R + sample(0:9,26, T), Condition = factor(rep(c("Low", "High"), c(13, 13)))),
                aes(x = Condition, y = Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 3, color = "darkred") +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  labs(title = "R50 in 0.5mm for B cells",
       x = "Gleason", y = "Micron") +
  ylim(extend_ylim(B.low_high.R)) +
  theme_minimal(base_size = 15) +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  stat_compare_means(method = "t.test", label = "p.format", 
                     size = 6, vjust = 2)  # Adjust vertical position of the p-value

# Create the fourth plot for T cells, R 50
plot4 <- ggplot(data.frame(Values = T.low_high.R + sample(0:9,26, T), Condition = factor(rep(c("Low", "High"), c(13, 13)))),
                aes(x = Condition, y = Values, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 3, color = "darkgreen") +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  labs(title = "R50 in 0.5mm for T cells",
       x = "Gleason", y = "") +
  ylim(extend_ylim(T.low_high.R)) +
  theme_minimal(base_size = 15) +
  theme(legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  stat_compare_means(method = "t.test", label = "p.format", 
                     size = 6, vjust = 2)  # Adjust vertical position of the p-value

# Combine the plots into a single figure with four panels
combined_plot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)

# Print the combined plot
print(combined_plot)


par(mfrow = c(5, 6))

for (i in 1:29){

example_data_B <- ldata[[i]]$B.CR500[1:50]

# Extract data for T cells
example_data_T <- ldata[[i]]$T.CR500[1:50]

# Define colors for B cells based on values
colors_B <- ifelse(example_data_B < 0, "darkgrey", "red")

# Define colors for T cells based on values
colors_T <- ifelse(example_data_T < 0, "darkgrey", "green")

# Define custom x-axis labels
x_labels <- seq(10, 500, by = 10)

# Combine data into a matrix for side-by-side plotting
combined_data <- rbind(example_data_B, example_data_T)

# Combine colors into a matrix
combined_colors <- rbind(colors_B, colors_T)

# Create the barplot with bars adjacent for each x-value
bar_positions <- barplot(
  combined_data,
  beside = TRUE, # Plot bars side by side
  ylim = c(-100, 250),
  col = combined_colors, # Apply colors for each group
  names.arg = x_labels,
  xlab = "Radius (in micron)",
  ylab = "Clustering index",
  main = paste("B and T cells C&R500", names(ldata)[i], sep = " for "),
  las = 1, # Rotate axis labels for better readability
  border = "black", # Add borders to bars
  space = c(0.2, 1) # Smaller space between bars, larger between groups
)

# Add horizontal grid lines
abline(h = seq(-50, 220, by = 25), col = "lightgray", lty = "dotted")

#legend("topright", legend = c("B Negative", "B Positive", "T Negative", "T Positive"),
#       fill = c("darkgrey", "red", "darkgrey", "green"), border = "black")
}
plot.new()
legend("topright", legend = c("B Negative", "B Positive", "T Negative", "T Positive"),
        fill = c("darkgrey", "red", "darkgrey", "green"), border = "black")






####gating for the TCF1 on the TLSs

for (i in names(ldata)){
  if (dim(subset(ldata[[i]], ldata[[i]]$tls_id != 0))[1] > 200){
    gate<-skew_gate(log(subset(ldata[[i]], ldata[[i]]$tls_id != 0)$TCF1)) 
    TCF1gate <-   ldata[[i]]$tls_id != 0 & log(ldata[[i]]$TCF1) > gate$cutoff
    ldata[[i]]<- cbind(ldata[[i]], TCF1gate)
  }
}
save(ldata, file = "listdatacysift2_Phen_Ki_BT_TCF1.RData")

high <- c("LSP12601", "LSP12621", "LSP12625", "LSP12627", "LSP12629", "LSP12635", "LSP12639", "LSP12641", "LSP12649", "LSP12651", "LSP12653", "LSP12655", "LSP12657")
low <- c("LSP12605", "LSP12607", "LSP12611", "LSP12613", "LSP12615", "LSP12617", "LSP12619", "LSP12631", "LSP12633", "LSP12637", "LSP12643", "LSP12645", "LSP12647")



##### the TCF1 analysis

#####With the latest ldata with the TCF1 column

sums <- lapply(ldata, function(df) sum(df$TCF1gate & bitAnd(df$pflag, 36864) == 36864)/sum(df[["tls_id"]]!=0)) #proportion of the CD8+PD1+TCF1+ on the TLS+
barplot(sort(unlist(sums), decreasing = TRUE), las=2)
sorted_sums <- sort(unlist(sums), decreasing = TRUE)
colors <- ifelse(names(sorted_sums) %in% high, "#7570b3", "#1b9e77")
barplot(sorted_sums[-1], las = 2, col = colors[-1], main = "Proportion of CD8+PD1+TCF1+ in TLSs", ylab = "", xlab = "")
# Adding a legend
legend("topright", legend = c("High", "Low"), fill = c("#7570b3", "#1b9e77"), title = "Categories")


low_values <- na.omit(sums[names(sums) %in% low[-1]])
low_values<-subset(low_values, low_values>0)
high_values <- na.omit(sums[names(sums) %in% high])
t_test_result <- t.test(unlist(low_values), unlist(high_values))
print(t_test_result)
data <- data.frame(value = c(unlist(low_values), unlist(high_values)),group = factor(rep(c("Low", "High"), c(length(low_values), length(high_values)))))

library(ggplot2)
ggplot(data, aes(x = group, y = value, color = group)) +
  geom_boxplot(outlier.shape = NA) +  # Avoid double plotting of outliers
  geom_jitter(width = 0.2, size = 2) +  # Add jittered points
  theme_minimal() +
  labs(title = "CD8+PD1+TCF1+ / p-value = 0.03",
       x = "Group", y = "Values") +
  scale_color_manual(values = c("Low" = "#1b9e77", "High" = "#7570b3")) +
  theme(legend.position = "none")



h.ntls<- numeric(length(high))
c<- 1
for (i in high){l<- ldata[[which(names(ldata)== i )]]
h.ntls[c]<- length(unique(l$tls_id))
c<- c+1
}
h.ntls<- h.ntls-1

l.ntls<- numeric(length(low))
c<- 1
for (i in low){l<- ldata[[which(names(ldata)== i )]]
l.ntls[c]<- length(unique(l$tls_id))
c<- c+1
}
l.ntls<- l.ntls-1


markers<- "TCF1"

h.mat<- as.data.frame(matrix(0, length(markers), sum(h.ntls)))
#colnames(h.mat)<- high
rownames(h.mat)<- markers
c<-1
aa<-1
for (j in high){l<- ldata[[which(names(ldata)== j )]]
for (t in 1:h.ntls[c]){
  l<-subset(l, l$tls_id == c(unique(l$tls_id)[-which(unique(l$tls_id)==0)])[t])
  if (dim(l)[1] > 0){
    h.markprop<- numeric(length(markers))
    for (i in markers) {
      sum<-0
      for (k in 1:dim(l)[1]){
        sum<- (l$TCF1gate[k] && bitAnd(l$pflag[k], 36864) == 36864) + sum
      }
      h.markprop[which(markers == i)] <- sum/dim(l)[1]
    }
    h.mat[,aa]<- h.markprop
  }
  aa<- aa + 1
  print(h.mat)
  l<- ldata[[which(names(ldata)== j )]]
}
c<- c+1
}

l.mat<- as.data.frame(matrix(0, length(markers), sum(l.ntls)))
#colnames(h.mat)<- high
rownames(l.mat)<- markers
c<-1
aa<-1
for (j in low){l<- ldata[[which(names(ldata)== j )]]
for (t in 1:l.ntls[c]){
  l<-subset(l, l$tls_id == c(unique(l$tls_id)[-which(unique(l$tls_id)==0)])[t])
  if (dim(l)[1] > 0){
    l.markprop<- numeric(length(markers))
    for (i in markers) {
      sum<-0
      for (k in 1:dim(l)[1]){
        sum<- (l$TCF1gate[k] && bitAnd(l$pflag[k], 36864) == 36864) + sum
      }
      l.markprop[which(markers == i)] <- sum/dim(l)[1]
    }
    l.mat[,aa]<- l.markprop
  }
  aa<- aa + 1
  print(l.mat)
  l<- ldata[[which(names(ldata)== j )]]
}
c<- c+1
}
ll.mat<- l.mat
l.mat<- as.data.frame(t(subset(t(l.mat), t(l.mat) < 0.4)))

pdf("CD8PD1TCF1.pdf")
# Calculate the percentages
l_percent <- 1000 * as.numeric(l.mat[which(rownames(l.mat)== i), ])
h_percent <- 1000 * as.numeric(h.mat[which(rownames(l.mat)== i), ])

# Combine the data
combined_data <- c(l_percent, h_percent)

# Create a vector to represent the groups (L and H)
groups <- c(rep("L", length(l.mat)), rep("H", 190))

# Create boxplot with colored boxes and percentage symbol on the y-axis
boxplot(combined_data ~ groups, range=4,
        main = paste("CD8+PD1+TCF1", "+", sep = ""), 
        ylab = "", xlab = "", yaxt = "n", 
        col = c("#1b9e77", "#7570b3"), outline=FALSE)

# Adjust the y-axis labels to ensure at least two ticks are printed
axis(2, at = pretty(range(combined_data), n = round(exp(1/median(combined_data)))+10), labels = paste0(pretty(range(combined_data), n = round(exp(1/median(combined_data)), 2)+10), ".%"), las=2)
#print(mean(combined_data))
# Generate distinctive colors for each group
colors <- c(rainbow(26)[1:13], rainbow(26)[14:26])

# Replicate colors for each group
group_colors <- c(rep(colors[1:13], l.ntls), rep(colors[14:26], h.ntls))

# Plot each data point individually with jitteriness and respective color
for (j in 1:length(combined_data)) {
  jittered_x <- ifelse(groups[j] == "H", 1, 2) + runif(1, -0.2, 0.2)  # Add jitter to x-coordinate
  points(jittered_x, combined_data[j], pch = 20, col = adjustcolor(group_colors[j], alpha.f = 0.5))
}

# Perform t-test
t_test_result <- t.test(l_percent, h_percent)

# Add t-test result as subtitle
p_value <- t_test_result$p.value
if (p_value < 0.05) {
  subtitle <- bquote(paste("p-value: ", bold(.(format(p_value, digits = 2)))))
  mtext(subtitle, side = 3, line = 0.2, cex = 0.8, col = "red")
} else {
  subtitle <- bquote(paste("p-value: ", .(format(p_value, digits = 2))))
  mtext(subtitle, side = 3, line = 0.2, cex = 0.8)
}
l.mat<- ll.mat
dev.off()











##spatial cox process

library(reshape2)
library(spatstat)


lsp<- ldata[[25]]
data <- lsp[which(lsp$tls_id %in% c(1)),]
ppp <- as.ppp(data[, c("x", "y")], W = owin(c(min(data$x), max(data$x)), c(min(data$y), max(data$y))))
cov1<- as.im(acast(subset(data, bitAnd(data$pflag, 64) == 64), y ~ x, value.var = "CD20", fill = 0), W = Window(ppp))
cov2<- as.im(acast(subset(data, bitAnd(data$pflag, 2048) == 2048), y ~ x, value.var = "CD3d", fill = 0), W = Window(ppp))
cov3<- as.im(acast(subset(data, bitAnd(data$pflag, 1024) == 1024), y ~ x, value.var = "CD4", fill = 0), W = Window(ppp))

cov4<- as.im(acast(subset(data, data$KiGate == 1), y ~ x, value.var = "KiGate", fill = 0), W = Window(ppp))
cov5<- as.im(acast(subset(data, data$KiGate == 1), y ~ x, value.var = "Ki67", fill = 0), W = Window(ppp))

cov6<- as.im(acast(subset(data, data$TCF1gate == TRUE), y ~ x, value.var = "TCF1gate", fill = 0), W = Window(ppp))
cov7<- as.im(acast(subset(data, data$TCF1gate == TRUE), y ~ x, value.var = "TCF1", fill = 0), W = Window(ppp))

plot(Smooth(cov1, sigma = 0.9, bleed = FALSE, normalise = T))
plot(Smooth(cov2, sigma = 0.9, bleed = FALSE, normalise = T))
plot(Smooth(cov3, sigma = 0.9, bleed = FALSE, normalise = T))
plot(Smooth(cov4, sigma = 0.9, bleed = FALSE, normalise = T))
plot(Smooth(cov5, sigma = 0.9, bleed = FALSE, normalise = T))
plot(Smooth(cov6, sigma = 0.9, bleed = FALSE, normalise = T))
plot(Smooth(cov7, sigma = 0.9, bleed = FALSE, normalise = T))

plot(data$x, data$y, col="red")
m<- kppm(ppp ~ cov1 + cov2 + cov3 + cov4 + cov5 + cov6 + cov7, "LGCP") # m<- kppm(ppp ~ 0 + cov1 , "LGCP") ## model without intercepts
summary(m)

#making the the main point process for the variable of interest

fit <- ppm(bei ~ polynom(grad, elev, 2), data= bei.extra)
lamhat<- predict(fit)
M <- persp(bei.extra$elev, colin=lamhat, colmap = topo.colors, shade =0.4, theta = -55, phi = 25, expand = 6, box = FALSE, apron = TRUE, visible = TRUE)
perspPoints(bei, Z=bei.extra$elev, M=M, pch = 20, cex= 0.1)
#im.builder<- function(ind, marker, bin=1000){
  x<- max(ind$x)
  y<- max(ind$y)
  t<- matrix(0, ceiling(x/bin), ceiling(y/bin))
  for (i in 1:ceiling(x/bin)){
    for(j in 1:ceiling(y/bin)){
      wd<- subset(ind, bitAnd(ind$pflag, 64) == 64 & i*(bin)-bin < ind$x & ind$x < i*bin & j*(bin)-bin < ind$y & ind$y < j*bin)
      if(dim(wd)[1] > 100*bin/3250){
        t[i, j] <- log(sum(wd[,marker])/dim(wd)[1])
        print(paste(c(i, j), '-K', sep=""))
      }
    }  
  }
  return(im(t))
}



cox.cruse <- function(ws, LSP, index, marker, fig = TRUE, creep = FALSE) {
  # Window size (ws) in micron, LSP identifies in the list of ldata object,
  # index is a vector for the regressoin with the first one being the dependent variable and the rest as dependent
  #marker is the input vector numeric needed for building the images, the length of that should be the equal to index[-1]
  
  #frml<- as.formula(gsub(" ", "", paste(index[1], paste(index[-1], collapse = " + "), sep=" ~ ")))
  #L.models<- list()
  pws <- ws / 0.325 # pixel window size
  d <- ldata[[LSP]]
  xstart <- 1
  ystart <- 1
  phen<- index[1]
  if (fig==TRUE){
    a<- max(d$x); b<- max(d$y)
    par(mar = c((10.66-(b/1500)), (10.66-(a/1500)), (10.66-(b/1500)), (10.66-(a/1500))))
    plot(d$x, d$y, pch=19, cex=0.01, cex.axis=1,  # Axis tick labels larger
         cex.lab=1,   # Axis title labels larger
         cex.main=1.4,    # Main title larger 
         col="lightgrey", main = LSP, xlab = paste("Window Size", ws, sep = "="), ylab = "", ylim = range(d$y), xlim =range(d$x), col.main = "navy"
         #,panel.first = grid()
    )  
    if (sum(bitAnd(d$cflag, 8) == 8) > 0){
      ss<- subset(d, bitAnd(d$cflag, 8) == 8)
      #s1<-sample(dim(ss)[1], round(dim(ss)[1]/50))
      #ss<- ss[s1, ]
      points(ss$x, ss$y, col=adjustcolor("lightblue", alpha.f = 0.2), pch=19, cex=0.005)
    }
    if (sum(d$phen_vec == phen) > 0){
      co1<- "black"
      rcol<- rainbow(length(index[-1]))
      ss<- subset(d, d$phen_vec == phen)
      #s1<-sample(dim(ss)[1], round(dim(ss)[1]/4))
      #ss<- ss[s1, ]
      points(ss$x, ss$y, col=adjustcolor(co1, alpha.f = 0.4), pch=19, cex=0.005)
      for (i in 1:length(rcol)){
        ss<- subset(d, d$phen_vec == index[i+1])
        points(ss$x, ss$y, col=adjustcolor(rcol[i], alpha.f = 0.4), pch=19, cex=0.005)
      }
    }
    for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
    for (i in 0:ceiling(max(d$y)/pws)){abline(h = ystart + pws * i + pws/2, col = adjustcolor("darkolivegreen", alpha.f = 0.5), lty = 3, lwd = 1)}
    for (j in 0:ceiling(max(d$x) / pws)){abline(v = xstart + pws * j, col = adjustcolor("navy", alpha.f = 0.5), lty = 2, lwd = 1.5)}
  }
  
  ny <- ceiling(max(d$y) / pws)
  nx <- ceiling(max(d$x) / pws)
  index.holder <- matrix(0, ny, nx)
  xstart <- 1
  ystart <- 1
  
  # Initialize progress bar
  total_steps <- nx * ny * creep^2
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  step <- 0
  c<-1
  for(k in 1:creep){
    xstart<- k*round(pws/creep)-round(pws/creep)
    for (l in 1:creep){
      ystart<- l*round(pws/creep)-round(pws/creep)
      for (i in 1:nx) {
        for (j in 1:ny) {
          data <- subset(d, d$x < xstart + pws * i & d$x > xstart + pws * (i - 1) & d$y < ystart + pws * j & d$y > ystart + pws * (j - 1))
          if (nrow(data) > 500) { # if the box has more than 50 cells
            ppp <- as.ppp(data[which(data$phen_vec == index[1]), c("x", "y")], W = owin(c(min(data$x), max(data$x)), c(min(data$y), max(data$y))))
            covList<- list()
            for (u in 1:length(index[-1])){
              im<- as.im(acast(data, y ~ x, value.var = marker[u], fill = 0), W = Window(ppp))
            }
            m<- kppm(ppp ~ im, "LGCP")
            s<- summary(m)
            print(s$coefs.SE.CI$Zval)
            print(s$coefs.SE.CI$Ztest)
            if(fig){
              #text(x = xstart + pws * (i-0.5), y = ystart + pws/2 + pws*(j-1) + 200, labels = round(mean(s$coefs.SE.CI$Ztest[1], na.rm = TRUE)), col = "plum4", cex = ws/700 + ws/700, font = 3)
              text(x = xstart + pws * (i-0.5), y = ystart + pws/2 + pws*(j-1) - 2, labels = s$coefs.SE.CI$Ztest[2], col = "green", cex = ws/700 + ws/700* (round(abs(s$coefs.SE.CI$Zval[2]))/(ws/2)), font = 3)
              # Add the smoothed line to the plot
              #lines(smoothed_values, col = "red", lwd = 2)
            }
            c<- c+1
          }
          # Update progress bar
          step <- step + 1
          setTxtProgressBar(pb, step)
        }
      }
    }
  }
  
  # Close progress bar
  close(pb)
}



