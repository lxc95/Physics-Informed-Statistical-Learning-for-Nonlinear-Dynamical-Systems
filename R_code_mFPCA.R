###############################################################POD30
start_time <- Sys.time()
library(fdapace)
library(R.matlab)

path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path(path, "all_condition_filed210914_POD30.mat")

data <- readMat(pathname)
data_set<-data[["force.allcondition"]]

tt = seq(1, 12000, by=1)
M = 12000
N = 35

L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(tt,N), t(data_set))

FPCAdense2 <- FPCA(L3$Ly, L3$Lt, optns = list(useBinnedData='OFF'))

loworder_force_phi = FPCAdense2$phi

filename <- paste("loworder_force_phi_POD30", ".mat", sep = "")
writeMat(filename, forceFPCAbasis = loworder_force_phi)

filename <- paste("loworder_force_mean_POD30", ".mat", sep = "")
writeMat(filename, forcemean = FPCAdense2$mu)

FPCAdense2$cumFVE
plot(FPCAdense2)

#############################################################################

library(fdapace)
library(R.matlab)

path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path(path, "all_condition_filed210914_POD30.mat")

data <- readMat(pathname)
data_set<-data[["displacement.allcondition"]]

tt = seq(1, 12000, by=1)
M = 12000
N = 35

L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(tt,N), t(data_set))

FPCAdense2 <- FPCA(L3$Ly, L3$Lt, optns = list(useBinnedData='OFF'))

loworder_displacement_phi = FPCAdense2$phi

filename <- paste("loworder_displacement_phi_POD30", ".mat", sep = "")
writeMat(filename, displacementFPCAbasis = loworder_displacement_phi)

filename <- paste("loworder_displacement_mean_POD30", ".mat", sep = "")
writeMat(filename, displacementmean = FPCAdense2$mu)

FPCAdense2$cumFVE
plot(FPCAdense2)
end_time <- Sys.time()
start_time-end_time