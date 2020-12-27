library(Matrix)
library(readr)
library(matlib)
library(plm)
library(stargazer)
library(pheatmap)
library(pracma)
setwd("~/Documents/DSE/Graph Theory/Assign_2")
data6 <- read_csv("data6.csv", col_types = cols(X1 = col_skip()))

rankMatrix(data6)

svd_data <- svd(data6)

# Plot heatmap of dataset
pheatmap(data6, how_colnames=T, show_colnames=T, cluster_rows = F, cluster_cols = T, display_numbers=T, legend=F, treeheight_col = 0)

# Compute SVD 
U <- svd_data$u
V <- svd_data$v
Vt <- t(svd_data$v)
D <- diag(svd_data$d, nrow = 5,ncol = 5)
d <- svd_data$d

pheatmap(D,cluster_rows = F,cluster_cols = F,legend = F,display_numbers = T)

# Compute SVD 
svd_data <- svd(data6, nu=2,nv=2)
U <- svd_data$u
V <- svd_data$v
Vt <- t(svd_data$v)
D <- diag(svd_data$d, nrow = 2,ncol = 2)

U <- as.matrix(U)

# Print Matrix in Latex 
stargazer(U,type = "latex")
stargazer(Vt,type = "latex")
stargazer(D,type = "latex")
stargazer(Mb, type = "latex")

Mb <- U%*%D%*%Vt
MbrankMatrix(Mb)

#Frobenius norm 
sqrt(sum((data6 - Mb)**2))

sqrt(sum(d**2))
sqrt(sum(d[1:2]**2))


library(ggplot2)

len_d<- c(1,2,3,4,5)
# We can clearly see how in the first two sigmas is concentrated the power of information 
ggplot()+
        geom_smooth(col="darkred")+
        theme_minimal()+
        geom_point(col= "darkred")+
        aes(log10(d),len_d, label = ("dfs"))+
        labs(x = "log of sigmas", y ='j') 

# the cumulative sum of sigmas over its actual total sum of sigmas. we can clearly
# see how just the first two sigmas cover more than 80% of the total information.
ggplot()+
        geom_smooth(col="darkred")+
        theme_minimal()+
        geom_point(col= "darkred")+
        aes(cumsum(d)/sum(d),len_d)+
        labs(x = "Cumulative sum of sigmas over the sum of sigmas", y="j",title = '') 
        

# The process of Singular Value Decomposition clearly indicates that exists only two main
# "concepts" or groups of customers, the first group is mainly concerned with what we can call Services
# i.e. the services speed and the availability of parking. The second group is mainly concerned with 
# the variety of courses. 
Mb<- `row.names<-`(Mb,c("1","2","3","4","5","6","7","8","9","10",'11','12','13','14','15'))


pheatmap(Mb,how_colnames=T, cluster_rows = T, cluster_cols = F, display_numbers=T, legend=F, show_rownames = T,how_rownames = T)

#let's look for the similarity of the other users:
Scores_users<-Mb%*%svd_data$v
Scores_users<- `row.names<-`(Scores_users,c("1","2","3","4","5","6","7","8","9","10",'11','12','13','14','15'))

pheatmap(Scores_users,cluster_rows = T,cluster_cols = F,display_numbers=T, legend=F, show_rownames = T,how_colnames=T, treeheight_row = 0 )
S <- Scores_users
# compute cos(theta) all users
# Compute and plot similarity matrix
similarity <- matrix(NA, nrow = 15, ncol = 15)
for (i in c(1:15)){
        for (j in c(1:15)){
                similarity[i, j]<-dot(S[i,],S[j,])/(sqrt(sum(S[i,]^2)*sum((S[j,])^2)))
        }
}


